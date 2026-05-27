#!/usr/bin/env python3
"""
Diagnostic plots ported from TSSr 0.99.6's ggplot2-based functions.
Uses plotnine (a faithful Python port of ggplot2 grammar of graphics) so
the call shape and visual output mirror the R originals.

Commands:
  tsspy plot pca   -- principal-component-analysis biplot of raw samples
                      (mirrors TSSr::plotTssPCA)
  tsspy plot iqw   -- per-sample interquantile-width histograms of clusters
                      (mirrors TSSr::plotInterQuantile)
  tsspy plot shape -- per-sample shape-score histograms (PSS / SI)
                      (mirrors TSSr::plotShape)
"""

from __future__ import annotations
import typer
import pandas as pd
import numpy as np
from pathlib import Path
from typing import List, Optional
import logging

from plotnine import (
    ggplot, aes, geom_point, geom_histogram, facet_wrap,
    theme_minimal, theme, element_text, labs, scale_fill_brewer, scale_color_brewer,
    coord_cartesian,
)

logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

__version__ = "0.21.0"
app = typer.Typer(help=f"Diagnostic plots (plotnine / ggplot2-style) (v{__version__})")

KEY_COLS = ['chr', 'pos', 'strand']


# -----------------------------------------------------------------------------
# plot pca  (TSSr::plotTssPCA)
# -----------------------------------------------------------------------------

@app.command("pca")
def pca_command(
    tss_table: Path = typer.Option(..., "-i", "--input",
                                   help="Raw (or processed) TSS table (chr/pos/strand/<samples>)"),
    output: Path = typer.Option("PCA_plot.pdf", "-o", "--output",
                                help="Output PDF/PNG (format inferred from extension)"),
    tss_threshold: int = typer.Option(10, "--tss-threshold",
                                      help="Keep TSS rows where any sample >= this value (TSSr default: 10)"),
    merge_labels: Optional[str] = typer.Option(None, "--merge-labels",
                                               help="Per-merged-group label (space-separated), e.g. 'YPD Arrest'"),
    merge_index: Optional[str] = typer.Option(None, "--merge-index",
                                              help="1-indexed group per sample (space-separated), e.g. '1 1 2 2'"),
    width: float = typer.Option(8.0, "--width"),
    height: float = typer.Option(8.0, "--height"),
):
    """
    PCA biplot of samples — TSSr::plotTssPCA.

    With default TSSr defaults the input is `@TSSrawMatrix` (pre-merge raw
    counts). Each row = TSS position, each sample column = read count;
    rows are filtered to those with any sample >= --tss-threshold, then
    PCA is run on the sample-by-position matrix (samples = observations).
    """
    from sklearn.decomposition import PCA

    df = pd.read_csv(tss_table, sep="\t")
    sample_cols = [c for c in df.columns if c not in KEY_COLS]
    if len(sample_cols) < 2:
        raise typer.BadParameter("PCA requires >= 2 sample columns.")

    sub = df[sample_cols]
    mask = (sub >= tss_threshold).any(axis=1)
    sub = sub[mask]
    if sub.empty:
        raise typer.BadParameter(f"No TSS rows after --tss-threshold {tss_threshold} filter.")
    logger.info(f"PCA on {sub.shape[1]} samples x {sub.shape[0]} TSS positions "
                f"(after threshold filter)")

    # PCA expects observations as rows. R does `t(tss)` (sample-by-position).
    X = sub.T.to_numpy(dtype=np.float64)
    pca = PCA(n_components=2)
    coords = pca.fit_transform(X)
    var_pct = pca.explained_variance_ratio_ * 100

    # Build label vector (per sample). R: `s <- sampleLabelsMerged[mergeIndex]`
    if merge_labels and merge_index:
        labels = merge_labels.split()
        idx = [int(x) for x in merge_index.split()]
        if len(idx) != len(sample_cols):
            raise typer.BadParameter(
                f"merge-index has {len(idx)} entries, sample columns are {len(sample_cols)}")
        label_per_sample = [labels[i - 1] for i in idx]
    else:
        # Default: each sample is its own group
        label_per_sample = list(sample_cols)

    plot_df = pd.DataFrame({
        "PC1": coords[:, 0],
        "PC2": coords[:, 1],
        "sample": sample_cols,
        "group": label_per_sample,
    })

    p = (ggplot(plot_df, aes("PC1", "PC2", color="group"))
         + geom_point(size=4)
         + labs(x=f"PC1 ({var_pct[0]:.1f}%)",
                y=f"PC2 ({var_pct[1]:.1f}%)",
                title="TSS PCA")
         + theme_minimal()
         + theme(text=element_text(size=12)))
    try:
        p = p + scale_color_brewer(type="qual", palette="Set1")
    except Exception:
        pass  # palette unavailable for some plotnine versions; defaults still fine

    p.save(str(output), width=width, height=height, verbose=False)
    typer.echo(f"PCA plot saved to {output}")


# -----------------------------------------------------------------------------
# plot iqw  (TSSr::plotInterQuantile)
# -----------------------------------------------------------------------------

@app.command("iqw")
def iqw_command(
    cluster_files: List[Path] = typer.Option(..., "-i", "--clusters",
                                             help="Per-sample cluster TSV (with columns tags + interquantile_width)"),
    sample_names: str = typer.Option(..., "-n", "--sample-names",
                                     help="Sample names, space-separated; one per --clusters file"),
    output: Path = typer.Option("Interquantile_plot.pdf", "-o", "--output"),
    tags_threshold: float = typer.Option(1.0, "--tags-threshold",
                                         help="Drop clusters with tags < threshold (TSSr default: 1)"),
    iqw_max: int = typer.Option(200, "--iqw-max",
                                help="Drop clusters with iqw > this (TSSr default: 200)"),
    bins: int = typer.Option(40, "--bins"),
    width: float = typer.Option(8.0, "--width"),
    height: float = typer.Option(6.0, "--height"),
):
    """Per-sample histograms of cluster interquantile-width (TSSr::plotInterQuantile)."""
    names = sample_names.split()
    if len(names) != len(cluster_files):
        raise typer.BadParameter(
            f"sample-names ({len(names)}) != cluster files ({len(cluster_files)})")

    frames = []
    for n, p in zip(names, cluster_files):
        df = pd.read_csv(p, sep="\t")
        df = df[(df["tags"] >= tags_threshold) & (df["interquantile_width"] <= iqw_max)]
        df = df[["interquantile_width"]].copy()
        df["sample"] = n
        frames.append(df)
    plot_df = pd.concat(frames, ignore_index=True)

    p = (ggplot(plot_df, aes("interquantile_width", fill="sample"))
         + geom_histogram(bins=bins, color="white", alpha=0.85)
         + facet_wrap("sample", scales="free_y")
         + labs(x="TC interquantile width q0.1-q0.9", y="Frequency",
                title="Interquantile width per sample")
         + theme_minimal()
         + theme(text=element_text(size=12), legend_position="none"))
    try:
        p = p + scale_fill_brewer(type="qual", palette="Set1")
    except Exception:
        pass

    p.save(str(output), width=width, height=height, verbose=False)
    typer.echo(f"IQW plot saved to {output}")


# -----------------------------------------------------------------------------
# plot shape  (TSSr::plotShape)
# -----------------------------------------------------------------------------

@app.command("shape")
def shape_command(
    shape_files: List[Path] = typer.Option(..., "-i", "--shape",
                                           help="Per-sample shape TSV (with column shape.score)"),
    sample_names: str = typer.Option(..., "-n", "--sample-names",
                                     help="Sample names, space-separated; one per --shape file"),
    output: Path = typer.Option("Shape_plot.pdf", "-o", "--output"),
    bins: int = typer.Option(40, "--bins"),
    width: float = typer.Option(8.0, "--width"),
    height: float = typer.Option(6.0, "--height"),
):
    """Per-sample histograms of cluster shape scores (TSSr::plotShape)."""
    names = sample_names.split()
    if len(names) != len(shape_files):
        raise typer.BadParameter(
            f"sample-names ({len(names)}) != shape files ({len(shape_files)})")

    frames = []
    for n, p in zip(names, shape_files):
        df = pd.read_csv(p, sep="\t")
        if "shape.score" not in df.columns:
            raise typer.BadParameter(f"{p} has no 'shape.score' column")
        sub = df[["shape.score"]].copy()
        sub["sample"] = n
        frames.append(sub)
    plot_df = pd.concat(frames, ignore_index=True)

    p = (ggplot(plot_df, aes("shape.score", fill="sample"))
         + geom_histogram(bins=bins, color="white", alpha=0.85)
         + facet_wrap("sample", scales="free_y")
         + labs(x="shape score", y="Frequency",
                title="Promoter shape score per sample")
         + theme_minimal()
         + theme(text=element_text(size=12), legend_position="none"))
    try:
        p = p + scale_fill_brewer(type="qual", palette="Set1")
    except Exception:
        pass

    p.save(str(output), width=width, height=height, verbose=False)
    typer.echo(f"Shape plot saved to {output}")


if __name__ == "__main__":
    app()
