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

from plotnine import (
    ggplot, aes, geom_point, geom_histogram, facet_wrap,
    theme_minimal, theme, element_text, labs, scale_fill_brewer, scale_color_brewer,
    coord_cartesian,
)

logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

__version__ = "0.24.0"
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
# plot tss  (TSSr::plotTSS) — coolbox-based
# -----------------------------------------------------------------------------

# Per-sample colour cycle (matches TSSr's rainbow() vibe)
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


def _write_sample_bigwig(tss_df: pd.DataFrame, sample: str, chrom_sizes: dict,
                         path: str, strand: str, negate: bool):
    """Write one BigWig per sample × strand. If `negate`, values become negative
    (used for minus strand so coolbox draws downward bars relative to baseline)."""
    import pyBigWig
    sub = tss_df[tss_df["strand"] == strand][["chr", "pos", sample]].copy()
    sub = sub[sub[sample] > 0]
    sub[sample] = sub[sample].astype(float)
    if negate:
        sub[sample] = -sub[sample]
    sub = sub[sub["chr"].isin(chrom_sizes)]
    present = [c for c in chrom_sizes if c in set(sub["chr"].unique())]
    chr_order = {c: i for i, c in enumerate(present)}
    sub["_o"] = sub["chr"].map(chr_order)
    sub = sub.sort_values(["_o", "pos"], kind="stable").reset_index(drop=True)

    bw = pyBigWig.open(path, "w")
    bw.addHeader([(c, chrom_sizes[c]) for c in present])
    if not sub.empty:
        for chrom, grp in sub.groupby("chr", sort=False):
            chroms = [chrom] * len(grp)
            starts = (grp["pos"].astype(int) - 1).tolist()
            ends = grp["pos"].astype(int).tolist()
            values = grp[sample].astype(float).tolist()
            bw.addEntries(chroms, starts, ends=ends, values=values)
    bw.close()


def _write_clusters_bed(cluster_df: pd.DataFrame, path: str):
    """Per-sample clusters → BED6 using q_0.1..q_0.9 as the span."""
    keep = cluster_df[["chr", "q_0.1", "q_0.9", "cluster", "strand"]].copy()
    keep["start_0based"] = (keep["q_0.1"].astype(int) - 1).clip(lower=0)
    keep["end"] = keep["q_0.9"].astype(int)
    keep["name"] = keep["cluster"].astype(int).astype(str)
    keep["score"] = 0
    keep[["chr", "start_0based", "end", "name", "score", "strand"]].to_csv(
        path, sep="\t", index=False, header=False)


def _write_gene_bed(ref_df: pd.DataFrame, path: str):
    """All gene rows → BED6 for the gene-model track."""
    out = ref_df[["chr", "start", "end", "gene_id", "strand"]].copy()
    out["start_0based"] = (out["start"].astype(int) - 1).clip(lower=0)
    out["end"] = out["end"].astype(int)
    out["score"] = 0
    out["name"] = out["gene_id"]
    out[["chr", "start_0based", "end", "name", "score", "strand"]].to_csv(
        path, sep="\t", index=False, header=False)


def _build_coolbox_frame(gene_bed_path: str,
                         sample_bw_plus_paths: List[str],
                         sample_bw_minus_paths: List[str],
                         sample_cluster_bed_paths: List[str],
                         sample_names: List[str]):
    """Compose coolbox Frame: XAxis + gene + (clusters + + strand bw + - strand bw) per sample."""
    from coolbox.api import XAxis, BED, BigWig

    frame = XAxis()
    frame = frame + BED(gene_bed_path, gene_style="flybase", title="gene",
                        height=1.0, color="#4682B4", labels=True, fontsize=8)
    for i, (sample, bw_plus, bw_minus, clu) in enumerate(zip(
            sample_names, sample_bw_plus_paths, sample_bw_minus_paths,
            sample_cluster_bed_paths)):
        color = _TRACK_PALETTE[i % len(_TRACK_PALETTE)]
        frame = (frame
                 + BED(clu, title=f"{sample} clusters", color=color,
                       height=0.6, labels=True, fontsize=7, gene_style="normal")
                 + BigWig(bw_plus, title=f"{sample} TSS (+)", color=color,
                          height=1.0, style="fill", line_width=0.6)
                 + BigWig(bw_minus, title=f"{sample} TSS (-)", color=color,
                          height=1.0, style="fill", line_width=0.6))
    return frame


def _plot_one_gene_coolbox(gene_row,
                           sample_names: List[str],
                           sample_bw_plus_paths: List[str],
                           sample_bw_minus_paths: List[str],
                           sample_cluster_bed_paths: List[str],
                           gene_bed_path: str,
                           up_dis: int, down_dis: int):
    """Render one gene's TSS plot via coolbox; return matplotlib Figure."""
    gene_chr = gene_row["chr"]
    gene_start = int(gene_row["start"])
    gene_end = int(gene_row["end"])
    gene_strand = gene_row["strand"]
    gene_id = gene_row.get("gene_id", "?")

    if gene_strand == "+":
        x_lo = max(1, gene_start - up_dis)
        x_hi = gene_end + down_dis
    else:
        x_lo = max(1, gene_start - down_dis)
        x_hi = gene_end + up_dis

    frame = _build_coolbox_frame(gene_bed_path, sample_bw_plus_paths,
                                  sample_bw_minus_paths,
                                  sample_cluster_bed_paths, sample_names)
    fig = frame.plot(f"{gene_chr}:{x_lo}-{x_hi}", close_fig=False)
    fig.suptitle(f"{gene_id}  ({gene_chr}:{gene_start}-{gene_end} {gene_strand})",
                 fontsize=12, y=1.02)
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
    reference: Optional[Path] = typer.Option(None, "-r", "--reference",
                                             help="Reference FASTA (chrom sizes for BigWig)"),
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
    import tempfile
    from TSSpy.bigwig import _read_chrom_sizes

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
    if reference is None:
        raise typer.BadParameter("--reference FASTA is required (for BigWig chrom sizes)")

    chrom_sizes = _read_chrom_sizes(str(reference), None)

    # Prepare temp files for coolbox: per-sample BigWig, per-sample cluster BED, gene BED
    with tempfile.TemporaryDirectory(prefix="tsspy_plottss_") as tmpdir:
        gene_bed_path = f"{tmpdir}/genes.bed"
        _write_gene_bed(ref, gene_bed_path)

        sample_bw_plus_paths = []
        sample_bw_minus_paths = []
        sample_cluster_bed_paths = []
        for sample, clu_path in zip(names, cluster_inputs):
            bw_p = f"{tmpdir}/{sample}.plus.bw"
            bw_m = f"{tmpdir}/{sample}.minus.bw"
            _write_sample_bigwig(tss_df, sample, chrom_sizes, bw_p, "+", negate=False)
            _write_sample_bigwig(tss_df, sample, chrom_sizes, bw_m, "-", negate=True)
            sample_bw_plus_paths.append(bw_p)
            sample_bw_minus_paths.append(bw_m)
            clu_df = pd.read_csv(clu_path, sep="\t")
            cb_path = f"{tmpdir}/{sample}.clusters.bed"
            _write_clusters_bed(clu_df, cb_path)
            sample_cluster_bed_paths.append(cb_path)

        with PdfPages(str(output)) as pdf:
            for _, row in ref_sub.iterrows():
                fig = _plot_one_gene_coolbox(
                    row, names, sample_bw_plus_paths, sample_bw_minus_paths,
                    sample_cluster_bed_paths, gene_bed_path,
                    up_dis, down_dis)
                pdf.savefig(fig, bbox_inches="tight")
                plt.close(fig)
    typer.echo(f"Wrote {len(ref_sub)}-page TSS browser PDF to {output}")


if __name__ == "__main__":
    app()
