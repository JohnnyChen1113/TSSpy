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

import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.patches import FancyArrow, Rectangle
from plotnine import (
    ggplot, aes, geom_point, geom_histogram, facet_wrap,
    theme_minimal, theme, element_text, labs, scale_fill_brewer, scale_color_brewer,
    coord_cartesian,
)

logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

__version__ = "0.22.0"
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


# -----------------------------------------------------------------------------
# plot tss  (TSSr::plotTSS)
# -----------------------------------------------------------------------------

# Per-sample bar/cluster track palette — cycled like R's rainbow()
_TRACK_PALETTE = ["#FF0000", "#0000FF", "#008000", "#FFA500",
                  "#800080", "#00CED1", "#A52A2A", "#FF1493"]


def _load_gene_ref(ref_table: Optional[Path], gff: Optional[Path]) -> pd.DataFrame:
    """Resolve gene reference table from either explicit TSV or GFF parser."""
    if ref_table is not None:
        df = pd.read_csv(ref_table, sep="\t")
        if "seqnames" in df.columns and "chr" not in df.columns:
            df = df.rename(columns={"seqnames": "chr"})
        return df
    if gff is not None:
        from TSSpy.gene_assign import load_genes_from_gff
        df = load_genes_from_gff(str(gff))
        df = df.rename(columns={"seqnames": "chr"})
        return df
    raise typer.BadParameter("Need --ref-table or --annotation (GFF).")


def _plot_one_gene_tss(gene_row, tss_window: pd.DataFrame,
                       per_sample_clusters: dict, sample_cols: List[str],
                       up_dis: int, down_dis: int, bidirection: bool, y_fixed: bool):
    """Render one gene's TSS browser-style figure to a matplotlib Figure."""
    gene_chr = gene_row["chr"]
    gene_start = int(gene_row["start"])
    gene_end = int(gene_row["end"])
    gene_strand = gene_row["strand"]
    gene_id = gene_row.get("gene_id", "?")

    if gene_strand == "+":
        x_lo = gene_start - up_dis
        x_hi = gene_end + down_dis
    else:
        x_lo = gene_start - down_dis
        x_hi = gene_end + up_dis

    n_samples = len(sample_cols)
    # Rows: 1 axis, 1 gene model, n_samples * 2 (cluster + TSS)
    n_rows = 2 + 2 * n_samples
    height_ratios = [0.6, 0.8] + [0.5, 1.5] * n_samples
    fig, axes = plt.subplots(n_rows, 1, figsize=(10, 1.5 + 0.6 * n_rows),
                             sharex=True,
                             gridspec_kw={"height_ratios": height_ratios})
    fig.suptitle(f"{gene_id}  ({gene_chr}:{gene_start}-{gene_end} {gene_strand})", fontsize=12)

    # Row 0: genomic axis
    ax_axis = axes[0]
    ax_axis.set_yticks([])
    ax_axis.set_frame_on(False)
    ax_axis.set_xlim(x_lo, x_hi)
    ax_axis.tick_params(axis="x", top=True, labeltop=True, bottom=False, labelbottom=False)
    ax_axis.text(x_lo, 0.5, f"{gene_chr}:{x_lo:,}-{x_hi:,}", fontsize=8,
                 va="center", ha="left", transform=ax_axis.transData)

    # Row 1: gene model (arrow)
    ax_gene = axes[1]
    ax_gene.set_yticks([])
    ax_gene.set_frame_on(False)
    ax_gene.set_xlim(x_lo, x_hi)
    arrow_y = 0.5
    arrow_height = 0.4
    if gene_strand == "+":
        ax_gene.add_patch(FancyArrow(
            gene_start, arrow_y, gene_end - gene_start, 0,
            width=arrow_height, length_includes_head=True,
            head_width=arrow_height * 1.5,
            head_length=min(150, (gene_end - gene_start) * 0.2),
            facecolor="#4682B4", edgecolor="black", linewidth=0.5))
    else:
        ax_gene.add_patch(FancyArrow(
            gene_end, arrow_y, gene_start - gene_end, 0,
            width=arrow_height, length_includes_head=True,
            head_width=arrow_height * 1.5,
            head_length=min(150, (gene_end - gene_start) * 0.2),
            facecolor="#4682B4", edgecolor="black", linewidth=0.5))
    ax_gene.text(0.0, 0.5, "gene", transform=ax_gene.transAxes,
                 ha="right", va="center", fontsize=8)
    ax_gene.set_ylim(0, 1)

    # Determine shared y range for TSS tracks if y_fixed
    if y_fixed:
        all_vals = tss_window[sample_cols].to_numpy()
        if all_vals.size:
            y_max = float(np.nanmax(np.abs(all_vals))) * 1.1
        else:
            y_max = 1.0
    else:
        y_max = None

    for idx, sample in enumerate(sample_cols):
        color = _TRACK_PALETTE[idx % len(_TRACK_PALETTE)]
        # Row 2 + 2*idx: cluster boxes
        ax_clu = axes[2 + 2 * idx]
        ax_clu.set_yticks([])
        ax_clu.set_xlim(x_lo, x_hi)
        ax_clu.set_frame_on(False)
        ax_clu.text(0.0, 0.5, f"{sample}\nclusters", transform=ax_clu.transAxes,
                    ha="right", va="center", fontsize=8)
        clu = per_sample_clusters.get(sample)
        if clu is not None:
            sub = clu[(clu["chr"] == gene_chr) &
                      (clu["strand"] == gene_strand) &
                      (clu["q_0.1"] >= x_lo) & (clu["q_0.9"] <= x_hi)]
            for _, r in sub.iterrows():
                ax_clu.add_patch(Rectangle((r["q_0.1"], 0.2), r["q_0.9"] - r["q_0.1"], 0.6,
                                           facecolor=color, alpha=0.4, edgecolor=color))
                if "cluster" in r and pd.notna(r["cluster"]):
                    ax_clu.text((r["q_0.1"] + r["q_0.9"]) / 2, 0.5, str(int(r["cluster"])),
                                ha="center", va="center", fontsize=7, color="black")
        ax_clu.set_ylim(0, 1)

        # Row 3 + 2*idx: TSS bars
        ax_tss = axes[3 + 2 * idx]
        ax_tss.set_xlim(x_lo, x_hi)
        ax_tss.spines["top"].set_visible(False)
        ax_tss.spines["right"].set_visible(False)
        ax_tss.text(0.0, 0.5, f"{sample}\nTSS (TPM)", transform=ax_tss.transAxes,
                    ha="right", va="center", fontsize=8)
        ax_tss.axhline(0, color="grey", linewidth=0.5)
        if bidirection:
            sub = tss_window[(tss_window["chr"] == gene_chr) &
                             (tss_window["pos"] >= x_lo) & (tss_window["pos"] <= x_hi)]
        else:
            sub = tss_window[(tss_window["chr"] == gene_chr) &
                             (tss_window["strand"] == gene_strand) &
                             (tss_window["pos"] >= x_lo) & (tss_window["pos"] <= x_hi)]
        positions = sub["pos"].to_numpy()
        values = sub[sample].astype(float).to_numpy()
        # Minus-strand TSS values are negated (TSSr convention; we mirror that)
        minus_mask = (sub["strand"].to_numpy() == "-")
        plot_vals = values.copy()
        plot_vals[minus_mask] = -np.abs(plot_vals[minus_mask])
        plot_vals[~minus_mask] = np.abs(plot_vals[~minus_mask])
        ax_tss.vlines(positions, 0, plot_vals, color=color, linewidth=0.8)
        if y_max is not None and y_max > 0:
            ax_tss.set_ylim(-y_max, y_max)

    axes[-1].set_xlabel(f"{gene_chr} position")
    plt.tight_layout(rect=(0, 0, 1, 0.96))
    return fig


@app.command("tss")
def tss_command(
    tss_table: Path = typer.Option(..., "-t", "--tss",
                                   help="TSS table (chr/pos/strand/<samples>, e.g. filterTSS output)"),
    cluster_inputs: List[Path] = typer.Option(..., "-c", "--clusters",
                                              help="Per-sample cluster file (consensus or tag)"),
    sample_names: str = typer.Option(..., "-n", "--sample-names",
                                     help="Sample names, space-separated; one per --clusters"),
    annotation: Optional[Path] = typer.Option(None, "-a", "--annotation",
                                              help="GFF3 annotation (one row per gene)"),
    ref_table: Optional[Path] = typer.Option(None, "--ref-table",
                                             help="Alternative to --annotation: pre-parsed gene TSV"),
    genes: str = typer.Option(..., "--genes",
                              help='Gene IDs to plot, space-separated (e.g. "YAL001C YPK1")'),
    output: Path = typer.Option("TSS_graphs.pdf", "-o", "--output",
                                help="Multi-page PDF output (one page per gene)"),
    up_dis: int = typer.Option(500, "--up-dis",
                               help="Bp upstream of gene to include (TSSr default: 500)"),
    down_dis: int = typer.Option(500, "--down-dis",
                                 help="Bp downstream of gene to include"),
    bidirection: bool = typer.Option(True, "--bidirection/--no-bidirection",
                                     help="Include TSSs on the opposite strand in the same window"),
    y_fixed: bool = typer.Option(True, "--y-fixed/--no-y-fixed",
                                 help="Use one shared y-range for all samples (TSSr default: TRUE)"),
):
    """
    Browser-style TSS plot per gene (TSSr::plotTSS).

    Multi-track stacked layout per page: genomic axis + gene model arrow
    + (cluster boxes + TSS bars) per sample.
    """
    names = sample_names.split()
    if len(names) != len(cluster_inputs):
        raise typer.BadParameter(
            f"sample-names ({len(names)}) != clusters ({len(cluster_inputs)})")
    requested_genes = genes.split()

    ref = _load_gene_ref(ref_table, annotation)
    ref_sub = ref[ref["gene_id"].isin(requested_genes)]
    missing = set(requested_genes) - set(ref_sub["gene_id"])
    if missing:
        logger.warning(f"requested genes not found in annotation: {sorted(missing)}")
    if ref_sub.empty:
        raise typer.BadParameter("None of the requested genes were found.")

    tss_df = pd.read_csv(tss_table, sep="\t")
    missing_cols = [n for n in names if n not in tss_df.columns]
    if missing_cols:
        raise typer.BadParameter(f"sample cols missing in TSS table: {missing_cols}")

    per_sample_clusters = {n: pd.read_csv(p, sep="\t") for n, p in zip(names, cluster_inputs)}

    with PdfPages(str(output)) as pdf:
        for _, row in ref_sub.iterrows():
            fig = _plot_one_gene_tss(row, tss_df, per_sample_clusters, names,
                                     up_dis, down_dis, bidirection, y_fixed)
            pdf.savefig(fig)
            plt.close(fig)
    typer.echo(f"Wrote {len(ref_sub)}-page TSS browser PDF to {output}")


if __name__ == "__main__":
    app()
