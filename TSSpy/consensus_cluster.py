#!/usr/bin/env python3
"""
Cross-sample consensus clusters — strict parity with TSSr 0.99.6 consensusCluster().

Mirrors:
  TSSr::consensusCluster()        (ConsensusMethods.R)
  TSSr::.getConsensus()           (ConsensusFunctions.R)
  TSSr::.getConsensusQuantile()   (ConsensusFunctions.R)

Algorithm:
  Phase 1 — build consensus GRanges set across all samples
    For each sample's tag clusters, construct fixed-width windows of
    [dominant_tss - round(dis/2), dominant_tss + round(dis/2)] (both inclusive).
    Sample 1: self-union (reduce) into disjoint ranges.
    Sample i (i >= 2): findOverlaps with current consensus, then concatenate
      union(overlap_pairs) + non-overlapping_from_1 + non-overlapping_from_2.
    Sort by (strand, chr, start) and assign 1-based consensusCluster ID.

  Phase 2 — per-sample quantile reconstruction
    For each consensus range gr[x]:
      Find this sample's tag clusters whose dominant_tss is in [gr.start, gr.end].
      If any: pull this sample's TSS positions in [min(tc.start), max(tc.end)],
              recompute tags-sum, dominant_tss, q_0.1, q_0.9, interquantile_width.
      If none: omit this sample's row for this consensus cluster.
    Sort per-sample output by (strand, chr, start).
"""

from __future__ import annotations
import typer
import pandas as pd
import numpy as np
from pathlib import Path
from typing import List, Optional, Dict, Tuple
import logging

logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

__version__ = "0.16.0"
app = typer.Typer(help=f"Cross-sample consensus clustering (v{__version__})")

KEY_COLS = ['chr', 'pos', 'strand']
OUTPUT_COLS = ['cluster', 'chr', 'start', 'end', 'strand',
               'dominant_tss', 'tags', 'tags.dominant_tss',
               'q_0.1', 'q_0.9', 'interquantile_width']


# -----------------------------------------------------------------------------
# GRanges-like primitives (per (chr, strand) stratification)
# -----------------------------------------------------------------------------

def _build_windows(tc: pd.DataFrame, dis: int) -> pd.DataFrame:
    """Construct ±round(dis/2) windows centered on each tc's dominant_tss."""
    half = int(round(dis / 2.0))  # banker's rounding matches R round()
    out = tc[['chr', 'strand', 'dominant_tss']].copy()
    out['start'] = out['dominant_tss'] - half
    out['end'] = out['dominant_tss'] + half
    return out[['chr', 'strand', 'start', 'end']]


def _grange_reduce_per_strata(ranges: pd.DataFrame) -> pd.DataFrame:
    """
    Per (chr, strand): sort by start, merge ranges with gap <= 0
    (overlapping or touching). Mirrors IRanges::reduce default behaviour
    (min.gapwidth = 1L means merge when gap < 1, i.e. gap == 0 too).
    """
    if ranges.empty:
        return ranges.copy()
    out = []
    for (c, s), grp in ranges.groupby(['chr', 'strand'], sort=False):
        grp_sorted = grp.sort_values('start').reset_index(drop=True)
        cur_start = int(grp_sorted.at[0, 'start'])
        cur_end = int(grp_sorted.at[0, 'end'])
        for i in range(1, len(grp_sorted)):
            s_i = int(grp_sorted.at[i, 'start'])
            e_i = int(grp_sorted.at[i, 'end'])
            if s_i <= cur_end + 1:  # overlap or touching => merge
                cur_end = max(cur_end, e_i)
            else:
                out.append((c, s, cur_start, cur_end))
                cur_start, cur_end = s_i, e_i
        out.append((c, s, cur_start, cur_end))
    return pd.DataFrame(out, columns=['chr', 'strand', 'start', 'end'])


def _find_overlaps(gr1: pd.DataFrame, gr2: pd.DataFrame) -> Tuple[np.ndarray, np.ndarray]:
    """
    Per (chr, strand): return (query_idx, subject_idx) arrays of overlapping pairs.
    Inclusive overlap: ranges [a,b] and [c,d] overlap iff max(a,c) <= min(b,d).
    """
    if gr1.empty or gr2.empty:
        return np.array([], dtype=np.int64), np.array([], dtype=np.int64)

    q_idx_list = []
    s_idx_list = []
    # Group both gr1 and gr2 by (chr, strand) to do O(N+M) sweep within strata.
    gr1 = gr1.reset_index(drop=False).rename(columns={'index': '_q'})
    gr2 = gr2.reset_index(drop=False).rename(columns={'index': '_s'})

    g1 = {k: v for k, v in gr1.groupby(['chr', 'strand'], sort=False)}
    g2 = {k: v for k, v in gr2.groupby(['chr', 'strand'], sort=False)}

    for key in g1:
        if key not in g2:
            continue
        a = g1[key].sort_values('start').reset_index(drop=True)
        b = g2[key].sort_values('start').reset_index(drop=True)
        # For each query in a, find subjects in b that overlap. Naive O(N*M) is fine
        # for the cluster-window sizes we deal with (thousands per strand max).
        a_starts = a['start'].values
        a_ends = a['end'].values
        b_starts = b['start'].values
        b_ends = b['end'].values
        a_qs = a['_q'].values
        b_ss = b['_s'].values
        for i in range(len(a)):
            s_a, e_a = a_starts[i], a_ends[i]
            # overlap: b.start <= e_a AND b.end >= s_a
            hits = np.where((b_starts <= e_a) & (b_ends >= s_a))[0]
            if len(hits) == 0:
                continue
            for h in hits:
                q_idx_list.append(int(a_qs[i]))
                s_idx_list.append(int(b_ss[h]))
    return np.array(q_idx_list, dtype=np.int64), np.array(s_idx_list, dtype=np.int64)


def _build_consensus_set(per_sample_tc: List[pd.DataFrame], dis: int) -> pd.DataFrame:
    """
    Build the cross-sample consensus GRanges set following TSSr's algorithm:
      gr <- union(gr1, gr1)
      for i in 2..N: gr <- .getConsensus(gr, cs[[i]], dis)
    The final set is NOT yet sorted or assigned IDs.
    """
    if not per_sample_tc:
        return pd.DataFrame(columns=['chr', 'strand', 'start', 'end'])

    # Sample 1: window construction + self-union (reduce)
    gr = _grange_reduce_per_strata(_build_windows(per_sample_tc[0], dis))

    for i in range(1, len(per_sample_tc)):
        gr2 = _build_windows(per_sample_tc[i], dis)
        q_idx, s_idx = _find_overlaps(gr, gr2)
        # union(gr[hit_q], gr2[hit_s])  -- with duplicates absorbed by reduce
        if len(q_idx) > 0:
            hit_q_unique = np.unique(q_idx)
            hit_s_unique = np.unique(s_idx)
            combined_hits = pd.concat([
                gr.iloc[hit_q_unique].assign(),
                gr2.iloc[hit_s_unique].assign(),
            ], ignore_index=True)
            unioned = _grange_reduce_per_strata(combined_hits)
        else:
            unioned = pd.DataFrame(columns=['chr', 'strand', 'start', 'end'])
            hit_q_unique = np.array([], dtype=np.int64)
            hit_s_unique = np.array([], dtype=np.int64)
        # gr[-hit_q] + gr2[-hit_s]
        non_hit_q = gr.drop(index=hit_q_unique) if len(hit_q_unique) else gr
        non_hit_s = gr2.drop(index=hit_s_unique) if len(hit_s_unique) else gr2
        gr = pd.concat([unioned, non_hit_q, non_hit_s], ignore_index=True)
    return gr


# -----------------------------------------------------------------------------
# Phase 2: per-sample quantile reconstruction
# -----------------------------------------------------------------------------

def _consensus_quantile_for_sample(consensus: pd.DataFrame,
                                   tc: pd.DataFrame,
                                   tss_sample: pd.DataFrame) -> pd.DataFrame:
    """
    For each consensus range, find sample's tcs whose dominant_tss lies in it.
    If any, pull TSS positions in [min(tc.start), max(tc.end)] and recompute stats.
    `tss_sample` is the sample's TSS table (chr, pos, strand, tags) filtered to tags > 0.
    """
    # Pre-index tc and tss for fast (chr, strand) lookup
    tc_by_cs = {k: v.sort_values('dominant_tss').reset_index(drop=True)
                for k, v in tc.groupby(['chr', 'strand'], sort=False)}
    tss_by_cs = {k: v.sort_values('pos').reset_index(drop=True)
                 for k, v in tss_sample.groupby(['chr', 'strand'], sort=False)}

    rows = []
    for cid, c_chr, c_strand, c_start, c_end in zip(
            consensus['consensusCluster'].values,
            consensus['chr'].values,
            consensus['strand'].values,
            consensus['start'].values,
            consensus['end'].values):
        key = (c_chr, c_strand)
        tc_grp = tc_by_cs.get(key)
        if tc_grp is None:
            continue
        # tc rows with dominant_tss in [c_start, c_end]
        mask = (tc_grp['dominant_tss'] >= c_start) & (tc_grp['dominant_tss'] <= c_end)
        temp = tc_grp[mask]
        if temp.empty:
            continue
        span_start = int(temp['start'].min())
        span_end = int(temp['end'].max())

        tss_grp = tss_by_cs.get(key)
        if tss_grp is None:
            continue
        s = tss_grp[(tss_grp['pos'] >= span_start) & (tss_grp['pos'] <= span_end)]
        if s.empty:
            continue
        s = s.sort_values('pos').reset_index(drop=True)

        # Integer-scaled tags (TPM rounded to 6 decimals -> scale by 1e6)
        tags_scaled = np.round(s['tags'].values * 1_000_000).astype(np.int64)
        total_scaled = int(tags_scaled.sum())
        tags_sum = total_scaled / 1_000_000.0

        # dominant_tss = pos of first max-tag row
        dom_idx = int(s['tags'].idxmax())
        dominant_tss = int(s.at[dom_idx, 'pos'])
        tags_dom = float(s.at[dom_idx, 'tags'])

        fwd = tags_scaled.cumsum()
        q1_idx = np.where(fwd * 10 > total_scaled)[0]
        q1 = int(s.at[int(q1_idx[0]), 'pos']) if len(q1_idx) else None

        rev = tags_scaled[::-1].cumsum()[::-1]
        q9_idx = np.where(rev * 10 > total_scaled)[0]
        q9 = int(s.at[int(q9_idx[-1]), 'pos']) if len(q9_idx) else None

        iqw = (q9 - q1 + 1) if (q1 is not None and q9 is not None) else 0

        rows.append({
            'cluster': int(cid),
            'chr': c_chr,
            'start': int(s['pos'].min()),
            'end': int(s['pos'].max()),
            'strand': c_strand,
            'dominant_tss': dominant_tss,
            'tags': tags_sum,
            'tags.dominant_tss': tags_dom,
            'q_0.1': q1,
            'q_0.9': q9,
            'interquantile_width': iqw,
        })

    if not rows:
        return pd.DataFrame(columns=OUTPUT_COLS)
    df = pd.DataFrame(rows)
    df = df.sort_values(['strand', 'chr', 'start']).reset_index(drop=True)
    return df[OUTPUT_COLS]


# -----------------------------------------------------------------------------
# Driver
# -----------------------------------------------------------------------------

def consensus_cluster(tss_df: pd.DataFrame,
                      per_sample_tc: Dict[str, pd.DataFrame],
                      dis: int = 50) -> Dict[str, pd.DataFrame]:
    """
    Run the full TSSr-style consensusCluster pipeline.
    Returns dict mapping sample name -> per-sample consensus DataFrame.
    """
    sample_order = list(per_sample_tc.keys())
    tc_list = [per_sample_tc[s] for s in sample_order]

    # Phase 1: cross-sample consensus
    gr = _build_consensus_set(tc_list, dis)
    if gr.empty:
        return {s: pd.DataFrame(columns=OUTPUT_COLS) for s in sample_order}

    gr = gr.sort_values(['strand', 'chr', 'start']).reset_index(drop=True)
    gr['consensusCluster'] = np.arange(1, len(gr) + 1, dtype=np.int64)

    # Phase 2: per-sample quantiles
    out: Dict[str, pd.DataFrame] = {}
    for s in sample_order:
        tc = per_sample_tc[s]
        tss_sample = tss_df[KEY_COLS + [s]].rename(columns={s: 'tags'}).copy()
        tss_sample = tss_sample[tss_sample['tags'] > 0]
        out[s] = _consensus_quantile_for_sample(gr, tc, tss_sample)
    return out


# -----------------------------------------------------------------------------
# CLI
# -----------------------------------------------------------------------------

@app.command("cluster")
def cluster_command(
    tss_input: Path = typer.Option(..., "-t", "--tss-input",
                                   help="Filtered TSS table from filterTSS"),
    cluster_inputs: List[Path] = typer.Option(..., "-i", "--input-clusters",
                                              help="Per-sample tag-cluster TSV files"),
    sample_names: str = typer.Option(..., "-n", "--sample-names",
                                     help='Sample names, space-separated; '
                                          'must match cluster file order AND a column in --tss-input'),
    output_prefix: Path = typer.Option(..., "-o", "--output-prefix",
                                       help="Output prefix; one file per sample"),
    dis: int = typer.Option(50, "-d", "--dis",
                            help="Window width around each dominant_tss; "
                                 "consensus uses ±round(dis/2). Default 50."),
):
    """Cross-sample consensus clusters (TSSr consensusCluster, dis=50 default)."""
    names = sample_names.split()
    if len(names) != len(cluster_inputs):
        raise typer.BadParameter(
            f"sample-names count ({len(names)}) != cluster files count ({len(cluster_inputs)})")

    tss_df = pd.read_csv(tss_input, sep='\t')
    missing = [n for n in names if n not in tss_df.columns]
    if missing:
        raise typer.BadParameter(f"sample columns not found in TSS table: {missing}")

    per_sample_tc: Dict[str, pd.DataFrame] = {}
    for n, p in zip(names, cluster_inputs):
        per_sample_tc[n] = pd.read_csv(p, sep='\t')
        logger.info(f"loaded {n}: {len(per_sample_tc[n])} tag clusters")

    out_per_sample = consensus_cluster(tss_df, per_sample_tc, dis=dis)
    for s, df in out_per_sample.items():
        out_path = Path(f"{output_prefix}.{s}.tsv")
        df.to_csv(out_path, sep='\t', index=False)
        logger.info(f"{s}: {len(df)} consensus rows -> {out_path}")


if __name__ == '__main__':
    app()
