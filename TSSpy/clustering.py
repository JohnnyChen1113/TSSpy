#!/usr/bin/env python3
"""
Cluster TSSs into tag clusters — strict parity with TSSr 0.99.6 clusterTSS().

Mirrors:
  TSSr::clusterTSS()      (ClusteringMethods.R)
  TSSr::.clusterByPeak()  (ClusteringFunctions.R)

Additionally supports `peakcluMax` (adamZhang_TSSr extension): greedy
non-maximum-suppression peak detection, re-validated by the standard
peakclu local-max test.

Input format: TSS table from filterTSS — chr, pos, strand, then one
column per merged sample (already TPM-normalized when fed via the
standard pipeline). Each sample column is processed independently and
written to a separate output file.
"""

from __future__ import annotations
import typer
import pandas as pd
import numpy as np
from pathlib import Path
from typing import List, Optional
import logging

logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

__version__ = "0.15.0"
app = typer.Typer(help=f"Cluster TSSs to infer core promoters (v{__version__})")

KEY_COLS = ['chr', 'pos', 'strand']
OUTPUT_COLS = ['cluster', 'chr', 'start', 'end', 'strand',
               'dominant_tss', 'tags', 'tags.dominant_tss',
               'q_0.1', 'q_0.9', 'interquantile_width']


# -----------------------------------------------------------------------------
# Core algorithm
# -----------------------------------------------------------------------------

def _rleid(values: np.ndarray) -> np.ndarray:
    """data.table-style run-length ID. values is a 2D array (n_rows × n_keys).
    A new run starts when ANY key changes from the previous row."""
    if len(values) == 0:
        return np.array([], dtype=np.int64)
    if values.ndim == 1:
        values = values.reshape(-1, 1)
    changes = np.any(values[1:] != values[:-1], axis=1)
    ids = np.empty(len(values), dtype=np.int64)
    ids[0] = 1
    ids[1:] = 1 + np.cumsum(changes)
    return ids


def _detect_peaks_peakclu(pos: np.ndarray, tags: np.ndarray, peak_distance: int) -> np.ndarray:
    """
    TSSr peakclu detection. For each TSS at index x (sorted by pos asc),
    peak[x] = x+1 if tags[x] equals max(tags) in the open window (pos-pd, pos+pd),
    where ties are broken by the first-encountered (lowest-pos) occurrence.
    Returns array of peak IDs (0 if not a peak; otherwise the 1-based row index).
    """
    n = len(pos)
    peak_id = np.zeros(n, dtype=np.int64)
    for x in range(n):
        p = pos[x]
        mask = (pos > p - peak_distance) & (pos < p + peak_distance)
        in_window_tags = tags[mask]
        in_window_pos = pos[mask]
        if in_window_tags.size == 0:
            continue
        max_tag = in_window_tags.max()
        first_max_pos = in_window_pos[in_window_tags == max_tag][0]
        if p == first_max_pos:
            peak_id[x] = x + 1
    return peak_id


def _detect_peaks_peakclumax(pos: np.ndarray, tags: np.ndarray, peak_distance: int) -> np.ndarray:
    """
    Adam Zhang's peakcluMax detection. Greedy: sort by tags descending,
    take highest, exclude [p-pd, p+pd] inclusive, repeat. Then re-validate
    each pick with the peakclu strict local-max test.
    Returns array of peak IDs.
    """
    n = len(pos)
    if n == 0:
        return np.zeros(0, dtype=np.int64)

    # Greedy max-first selection. Stable order: when tags tie, lower pos first.
    order = np.lexsort((pos, -tags))  # primary: -tags asc (= tags desc); secondary: pos asc
    available = np.ones(n, dtype=bool)
    peak_max_positions = []
    for idx in order:
        if not available[idx]:
            continue
        peak_max_positions.append(int(pos[idx]))
        # Remove [p-pd, p+pd] inclusive
        p = pos[idx]
        available &= ~((pos >= p - peak_distance) & (pos <= p + peak_distance))

    peak_max_set = set(peak_max_positions)

    # Re-validate with peakclu's strict local-max test
    peak_id = np.zeros(n, dtype=np.int64)
    for x in range(n):
        if int(pos[x]) not in peak_max_set:
            continue
        p = pos[x]
        mask = (pos > p - peak_distance) & (pos < p + peak_distance)
        in_window_tags = tags[mask]
        in_window_pos = pos[mask]
        if in_window_tags.size == 0:
            continue
        max_tag = in_window_tags.max()
        first_max_pos = in_window_pos[in_window_tags == max_tag][0]
        if p == first_max_pos:
            peak_id[x] = x + 1
    return peak_id


def _cluster_by_peak(group: pd.DataFrame,
                     peak_distance: int,
                     local_threshold: float,
                     extension_distance: int,
                     method: str) -> Optional[pd.DataFrame]:
    """
    Faithful port of TSSr's .clusterByPeak. Operates on a single (chr, strand)
    subset, already sorted by pos ascending, with columns: chr, pos, strand, tags.
    Returns a per-cluster DataFrame or None if empty.

    method ∈ {'peakclu', 'peakcluMax'}.
    """
    copied = group.copy().reset_index(drop=True)
    pos_arr = copied['pos'].to_numpy()
    tags_arr = copied['tags'].to_numpy()
    n = len(copied)
    if n == 0:
        return None

    # ---- peak detection ----
    if method == 'peakclu':
        peak_id = _detect_peaks_peakclu(pos_arr, tags_arr, peak_distance)
    elif method == 'peakcluMax':
        peak_id = _detect_peaks_peakclumax(pos_arr, tags_arr, peak_distance)
    else:
        raise ValueError(f"Unknown method {method!r}; expected 'peakclu' or 'peakcluMax'.")

    work = copied.copy()
    work['peak'] = peak_id
    work['ID'] = np.arange(1, n + 1)

    # ---- local filter (TSSr's asymmetric per-strand) ----
    # R's $tag partial-matches $tags (verified). The filter IS active.
    strand_val = group['strand'].iloc[0]
    drop_ids = set()
    for i in np.where(peak_id > 0)[0]:
        peak_pos = pos_arr[i]
        peak_tag = tags_arr[i]
        thresh = peak_tag * local_threshold
        if strand_val == '+':
            mask = (work['pos'] >= peak_pos) & (work['pos'] <= peak_pos + peak_distance)
        else:
            mask = (work['pos'] >= peak_pos - peak_distance) & (work['pos'] <= peak_pos)
        sub = work[mask]
        drop_ids.update(sub.loc[sub['tags'] < thresh, 'ID'].tolist())

    if drop_ids:
        work = work[~work['ID'].isin(drop_ids)].reset_index(drop=True)
        if len(work) == 0:
            return None
        pos_arr = work['pos'].to_numpy()
        tags_arr = work['tags'].to_numpy()
        peak_id = work['peak'].to_numpy()

    # ---- forward / reverse adjacency flags ----
    n = len(work)
    forward = np.zeros(n, dtype=np.int64)
    reverse = np.zeros(n, dtype=np.int64)
    if n > 1:
        forward[:-1] = (pos_arr[1:] < pos_arr[:-1] + extension_distance).astype(np.int64)
        reverse[1:] = (pos_arr[:-1] > pos_arr[1:] - extension_distance).astype(np.int64)
    work['forward'] = forward
    work['reverse'] = reverse

    # ---- rleid grouping on (peak, forward, reverse) ----
    rle = _rleid(work[['peak', 'forward', 'reverse']].to_numpy())
    work['rleid'] = rle
    collapsed = work.groupby('rleid', sort=True).agg(
        peak=('peak', 'max'),
        start=('pos', 'min'),
        end=('pos', 'max'),
        tags=('tags', 'sum'),
    ).reset_index()

    # ---- boundary extension (2-step lookahead each side, only for peak-containing groups) ----
    n_collapsed = len(collapsed)
    raw_bounds = []
    for x in range(n_collapsed):
        if collapsed.at[x, 'peak'] <= 0:
            continue
        start = int(collapsed.at[x, 'start'])
        end = int(collapsed.at[x, 'end'])
        # extend backward up to 2 non-peak groups
        if x - 1 >= 0 and collapsed.at[x - 1, 'peak'] <= 0 and \
                collapsed.at[x - 1, 'end'] > start - extension_distance:
            start = int(collapsed.at[x - 1, 'start'])
            if x - 2 >= 0 and collapsed.at[x - 2, 'peak'] <= 0 and \
                    collapsed.at[x - 2, 'end'] > start - extension_distance:
                start = int(collapsed.at[x - 2, 'start'])
        # extend forward up to 2 non-peak groups
        if x + 1 < n_collapsed - 1 and collapsed.at[x + 1, 'peak'] <= 0 and \
                collapsed.at[x + 1, 'start'] < end + extension_distance:
            end = int(collapsed.at[x + 1, 'end'])
            if x + 2 < n_collapsed - 1 and collapsed.at[x + 2, 'peak'] <= 0 and \
                    collapsed.at[x + 2, 'start'] < end + extension_distance:
                end = int(collapsed.at[x + 2, 'end'])
        raw_bounds.append([start, end])

    if not raw_bounds:
        return None
    bounds = pd.DataFrame(raw_bounds, columns=['V1', 'V2'])

    # ---- overlap merge: where V2[i] >= V1[i+1], merge ----
    # TSSr does this in a single pass: for each overlapping row, propagate V1[i] -> V1[i+1] then drop i.
    while True:
        overlap_idx = np.where(bounds['V2'].values[:-1] >= bounds['V1'].values[1:])[0]
        if len(overlap_idx) == 0:
            break
        # propagate the leftmost V1 to the next row, then drop the current row
        v1 = bounds['V1'].values.copy()
        for i in overlap_idx:
            v1[i + 1] = v1[i]
        bounds = bounds.assign(V1=v1).drop(index=overlap_idx).reset_index(drop=True)

    # ---- compute per-cluster metrics from the ORIGINAL un-filtered TSS data ----
    cluster_rows = []
    for i in range(len(bounds)):
        start = int(bounds.at[i, 'V1'])
        end = int(bounds.at[i, 'V2'])
        cluster_data = copied[(copied['pos'] >= start) & (copied['pos'] <= end)].copy()
        cluster_data = cluster_data.sort_values('pos').reset_index(drop=True)
        if cluster_data.empty:
            continue

        # dominant_tss = pos at first occurrence of max tags
        dominant_idx = int(cluster_data['tags'].idxmax())
        dominant_tss = int(cluster_data.at[dominant_idx, 'pos'])
        tags_dom = float(cluster_data.at[dominant_idx, 'tags'])

        # Quantiles + tags-sum via INTEGER arithmetic to avoid floating-point creep.
        # All tags are TSSr's round(x, 6) TPM values, so they're exact rationals
        # in 1e-6 units. Scale to integers, sum exactly, compare cumsum*10 vs total.
        tags_scaled = np.round(cluster_data['tags'].values * 1_000_000).astype(np.int64)
        total_scaled = int(tags_scaled.sum())
        tags_sum = total_scaled / 1_000_000.0

        # q_0.1: smallest pos where cumsum > 0.1 * total
        #         <=> 10 * cumsum > total
        fwd_cum_scaled = tags_scaled.cumsum()
        q1_idx = np.where(fwd_cum_scaled * 10 > total_scaled)[0]
        q1 = int(cluster_data.at[int(q1_idx[0]), 'pos']) if len(q1_idx) else None

        # q_0.9: max pos where reverse-cumsum > 0.1 * total
        rev_cum_scaled = tags_scaled[::-1].cumsum()[::-1]
        q9_idx = np.where(rev_cum_scaled * 10 > total_scaled)[0]
        q9 = int(cluster_data.at[int(q9_idx[-1]), 'pos']) if len(q9_idx) else None

        if q1 is None or q9 is None:
            iqw = 0
        else:
            iqw = q9 - q1 + 1

        cluster_rows.append({
            'chr': cluster_data['chr'].iloc[0],
            'start': start,
            'end': end,
            'strand': cluster_data['strand'].iloc[0],
            'dominant_tss': dominant_tss,
            'tags': tags_sum,
            'tags.dominant_tss': tags_dom,
            'q_0.1': q1,
            'q_0.9': q9,
            'interquantile_width': iqw,
        })

    if not cluster_rows:
        return None
    return pd.DataFrame(cluster_rows)


def cluster_one_sample(tss_df: pd.DataFrame,
                       sample_col: str,
                       peak_distance: int,
                       local_threshold: float,
                       extension_distance: int,
                       cluster_threshold: float,
                       method: str) -> pd.DataFrame:
    """Cluster a single sample column across all chrs and strands."""
    df = tss_df[KEY_COLS + [sample_col]].copy()
    df = df.rename(columns={sample_col: 'tags'})
    df = df[df['tags'] > 0]
    if df.empty:
        return pd.DataFrame(columns=OUTPUT_COLS)

    all_clusters = []
    for (chrom, strand), grp in df.groupby(['chr', 'strand'], sort=False):
        grp_sorted = grp.sort_values('pos').reset_index(drop=True)
        clusters = _cluster_by_peak(grp_sorted, peak_distance, local_threshold,
                                    extension_distance, method)
        if clusters is not None and not clusters.empty:
            all_clusters.append(clusters)

    if not all_clusters:
        return pd.DataFrame(columns=OUTPUT_COLS)

    out = pd.concat(all_clusters, ignore_index=True)
    # cluster threshold (TSSr does this BEFORE the final sort + .I assignment)
    out = out[out['tags'] > cluster_threshold].reset_index(drop=True)
    # sort by strand, chr, start (TSSr's setorder)
    out = out.sort_values(['strand', 'chr', 'start']).reset_index(drop=True)
    out['cluster'] = np.arange(1, len(out) + 1)
    return out[OUTPUT_COLS]


# -----------------------------------------------------------------------------
# CLI
# -----------------------------------------------------------------------------

@app.callback(invoke_without_command=True)
def main(
    ctx: typer.Context,
    input_file: Optional[Path] = typer.Option(None, "-i", "--input",
                                              help="Filtered TSS table (output of filterTSS)"),
    output_prefix: Optional[Path] = typer.Option(None, "-o", "--output-prefix",
                                                 help="Output prefix; per-sample files named <prefix>.<sample>.tsv"),
    samples: Optional[str] = typer.Option(None, "-s", "--samples",
                                          help="Comma-separated sample columns. Default: all non-key columns."),
    method: str = typer.Option("peakclu", "--method",
                               help="'peakclu' (TSSr default) or 'peakcluMax' (greedy NMS variant)"),
    peak_distance: int = typer.Option(100, "--peak-distance"),
    extension_distance: int = typer.Option(30, "--extension-distance"),
    local_threshold: float = typer.Option(0.02, "--local-threshold"),
    cluster_threshold: float = typer.Option(1.0, "--cluster-threshold"),
    version: bool = typer.Option(False, "--version", help="Show version and exit."),
):
    """Cluster TSSs into tag clusters (TSSr peakclu / peakcluMax)."""
    if version:
        typer.echo(f"clustering.py version {__version__}")
        raise typer.Exit(0)
    if input_file is None or output_prefix is None:
        typer.echo(ctx.get_help())
        raise typer.Exit(0)

    df = pd.read_csv(input_file, sep='\t')
    if samples:
        sample_cols = [s.strip() for s in samples.split(',')]
    else:
        sample_cols = [c for c in df.columns if c not in KEY_COLS]
    logger.info(f"Method: {method}  samples: {sample_cols}")
    logger.info(f"params: peakDistance={peak_distance}, extensionDistance={extension_distance}, "
                f"localThreshold={local_threshold}, clusterThreshold={cluster_threshold}")

    for sc in sample_cols:
        if sc not in df.columns:
            raise typer.BadParameter(f"sample column {sc!r} not in input")
        out = cluster_one_sample(df, sc, peak_distance, local_threshold,
                                 extension_distance, cluster_threshold, method)
        out_path = Path(f"{output_prefix}.{sc}.tsv")
        out.to_csv(out_path, sep='\t', index=False)
        logger.info(f"{sc}: {len(out):,} clusters -> {out_path}")


if __name__ == '__main__':
    app()
