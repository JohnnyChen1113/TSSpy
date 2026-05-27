#!/usr/bin/env python3
"""
Core-promoter shape scores — strict parity with TSSr 0.99.6 shapeCluster().

Mirrors TSSr::shapeCluster() (ShapeMethods.R).

For each cluster (from tagClusters or consensusClusters):
  Find this sample's TSS positions in [q_0.1, q_0.9] (inclusive).
  Let p_i = tags_i / sum(tags_i)
  PSS = -sum(p_i * log2(p_i)) * log2(interquantile_width)
  SI  = 2 + sum(p_i * log2(p_i))           (note: sum is negative)

No special-case branches: edge cases (singleton, iqw==1) fall out of the
formula naturally — log2(1) = 0 gives PSS=0 for iqw=1 or single-position
clusters, and gives SI=2 for single-position clusters. Matches TSSr.

The original cluster row is preserved verbatim; a `shape.score` column
(dot, matching TSSr) is appended. Row order is preserved (no resort).
"""

from __future__ import annotations
import math
import typer
import pandas as pd
import numpy as np
from pathlib import Path
from typing import List, Optional, Dict, Tuple
import logging

logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

__version__ = "0.17.0"
app = typer.Typer(help=f"Core-promoter shape scores PSS/SI (v{__version__})")

KEY_COLS = ['chr', 'pos', 'strand']


# -----------------------------------------------------------------------------
# Core algorithm
# -----------------------------------------------------------------------------

def _shape_score(tags: np.ndarray, iqw: int, method: str) -> float:
    """Per-cluster shape score, matching TSSr's formula directly.

    PSS = -sum_i (p_i * log2(p_i)) * log2(iqw),  p_i = tags_i / sum(tags)
    SI  = 2 + sum_i (p_i * log2(p_i))

    Edge cases fall out of the formula: singleton (single p=1) gives sum=0,
    so PSS=0 and SI=2. iqw=1 gives log2(1)=0, so PSS=0. Empty tags returns
    PSS=0, SI=2 (matches TSSr's empty-sapply behaviour where sum(numeric(0))=0).

    Floating-point note: R's log(x,2) is implemented as log(x)/log(2), which
    accumulates 1-2 ulps relative to IEEE log2(). We use R's formulation
    (np.log / LN2) for that reason. Even so, residual differences from R
    are bounded by ~1e-13 — well below the 6-decimal TPM input precision
    and any conceivable biological threshold.
    """
    if len(tags) == 0:
        return 0.0 if method == 'PSS' else 2.0
    total = float(tags.sum())
    if total <= 0:
        return 0.0 if method == 'PSS' else 2.0
    p = tags.astype(np.float64) / total
    LN2 = math.log(2.0)
    plogp_sum = float(np.sum(p * (np.log(p) / LN2)))
    if method == 'PSS':
        return -plogp_sum * (math.log(float(iqw)) / LN2)
    if method == 'SI':
        return 2.0 + plogp_sum
    raise ValueError(f"Unknown method {method!r}; expected 'PSS' or 'SI'.")


def add_shape_scores(cluster_df: pd.DataFrame,
                     tss_sample: pd.DataFrame,
                     method: str = 'PSS') -> pd.DataFrame:
    """
    Append a `shape.score` column to cluster_df.
    `tss_sample` must already be filtered to (tags > 0) for the relevant sample
    and have columns chr, pos, strand, tags.
    Preserves input row order (matches TSSr).
    """
    # Index TSS by (chr, strand) for fast lookup
    tss_by_cs = {k: v.sort_values('pos').reset_index(drop=True)
                 for k, v in tss_sample.groupby(['chr', 'strand'], sort=False)}

    chr_arr = cluster_df['chr'].to_numpy()
    strand_arr = cluster_df['strand'].to_numpy()
    q1_arr = cluster_df['q_0.1'].to_numpy()
    q9_arr = cluster_df['q_0.9'].to_numpy()
    iqw_arr = cluster_df['interquantile_width'].to_numpy()

    scores = np.empty(len(cluster_df), dtype=np.float64)
    for i in range(len(cluster_df)):
        sub = tss_by_cs.get((chr_arr[i], strand_arr[i]))
        if sub is None:
            tags = np.array([], dtype=np.float64)
        else:
            q1, q9 = q1_arr[i], q9_arr[i]
            mask = (sub['pos'] >= q1) & (sub['pos'] <= q9)
            tags = sub.loc[mask, 'tags'].to_numpy()
        scores[i] = _shape_score(tags, int(iqw_arr[i]), method)
    out = cluster_df.copy()
    out['shape.score'] = scores
    return out


def _extract_tss_for_sample(tss_df: pd.DataFrame, sample_col: str) -> pd.DataFrame:
    """Filter the full TSS table to one sample's positive-tag positions."""
    sub = tss_df[KEY_COLS + [sample_col]].copy()
    sub = sub.rename(columns={sample_col: 'tags'})
    sub = sub[sub['tags'] > 0]
    return sub


# -----------------------------------------------------------------------------
# CLI
# -----------------------------------------------------------------------------

@app.command("calculate")
def calculate_command(
    cluster_file: Path = typer.Option(..., "-c", "--clusters",
                                      help="Per-sample cluster file (consensus or tag)"),
    tss_file: Path = typer.Option(..., "-t", "--tss", help="Filtered TSS table"),
    output_file: Path = typer.Option(..., "-o", "--output"),
    sample: str = typer.Option(..., "-s", "--sample",
                               help="Sample column in --tss to use for this cluster file"),
    method: str = typer.Option("PSS", "-m", "--method", help="'PSS' or 'SI'"),
):
    """Shape score for one sample's clusters (TSSr shapeCluster, one sample)."""
    cluster_df = pd.read_csv(cluster_file, sep='\t')
    tss_df = pd.read_csv(tss_file, sep='\t')
    if sample not in tss_df.columns:
        available = [c for c in tss_df.columns if c not in KEY_COLS]
        raise typer.BadParameter(f"Sample '{sample}' not found. Available: {available}")
    tss_sample = _extract_tss_for_sample(tss_df, sample)
    out = add_shape_scores(cluster_df, tss_sample, method=method.upper())
    out.to_csv(output_file, sep='\t', index=False)
    logger.info(f"{sample}: {len(out)} clusters scored ({method}) -> {output_file}")


@app.command("batch")
def batch_command(
    cluster_inputs: List[Path] = typer.Option(..., "-c", "--clusters",
                                              help="Per-sample cluster TSV files"),
    sample_names: str = typer.Option(..., "-n", "--sample-names",
                                     help='Sample names, space-separated; must match TSS columns'),
    tss_file: Path = typer.Option(..., "-t", "--tss", help="Filtered TSS table"),
    output_prefix: Path = typer.Option(..., "-o", "--output-prefix",
                                       help="Output prefix; <prefix>.<sample>.tsv per sample"),
    method: str = typer.Option("PSS", "-m", "--method", help="'PSS' or 'SI'"),
):
    """Shape scores for all samples; mirrors TSSr shapeCluster's whole-object loop."""
    names = sample_names.split()
    if len(names) != len(cluster_inputs):
        raise typer.BadParameter(
            f"sample-names count ({len(names)}) != cluster files count ({len(cluster_inputs)})")
    tss_df = pd.read_csv(tss_file, sep='\t')
    missing = [n for n in names if n not in tss_df.columns]
    if missing:
        raise typer.BadParameter(f"sample columns missing from TSS table: {missing}")
    for n, p in zip(names, cluster_inputs):
        cluster_df = pd.read_csv(p, sep='\t')
        tss_sample = _extract_tss_for_sample(tss_df, n)
        out = add_shape_scores(cluster_df, tss_sample, method=method.upper())
        out_path = Path(f"{output_prefix}.{n}.tsv")
        out.to_csv(out_path, sep='\t', index=False)
        logger.info(f"{n}: {len(out)} clusters scored ({method}) -> {out_path}")


if __name__ == '__main__':
    app()
