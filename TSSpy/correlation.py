#!/usr/bin/env python3
"""
Sample correlation matrix + pairs-style plot — mirrors TSSr 0.99.6
plotCorrelation() and its internal `.plotCorrelation()` helper.

TSSr's reference behaviour (ExportFunctions.R:1-22):
  - Operates on @TSSrawMatrix (raw counts pre-normalization) by default.
  - `pairs(z, lower.panel=scatter, upper.panel=text_r, log="xy")`
      scatter: pch=".", col="#00AFBB"
      text:    r rounded to 2 dp, font size scaled by |r|
  - Diagonal: variable names (R `pairs()` default).
"""

from __future__ import annotations
import typer
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt


__version__ = "0.20.0"

# TSSr's scatter colour (`#00AFBB` from .plotCorrelation upper.panel)
TSSR_POINT_COLOR = "#00AFBB"
KEY_COLS = ["chr", "pos", "strand"]


def correlation(
    tss_table: str = typer.Option(..., "-i", "--input",
                                  help="TSS table (tab-delimited, columns: chr, pos, strand, <samples>)"),
    output: str = typer.Option(None, "-o", "--output",
                               help="Output correlation matrix CSV (optional)"),
    plot: bool = typer.Option(False, "--plot/--no-plot",
                              help="Also generate the R-pairs-style plot"),
    plot_file: str = typer.Option("correlation_pairs.png", "--plot-file",
                                  help="Output plot file (PNG / PDF; format inferred from extension)"),
    source: str = typer.Option("raw", "--source",
                               help="'raw' (TSSr-style log-log scatter) or 'processed' (log10(x+1) on linear axes)"),
    version: bool = typer.Option(False, "--version", help="Show version and exit"),
):
    """Pairwise sample correlation matrix + optional plot (TSSr-style)."""
    if version:
        typer.echo(f"TSSpy correlation v{__version__}")
        raise typer.Exit()

    df = pd.read_csv(tss_table, sep="\t")
    sample_cols = [c for c in df.columns if c not in KEY_COLS]
    if len(sample_cols) < 2:
        raise typer.BadParameter(f"need >=2 sample columns; got {sample_cols}")
    data = df[sample_cols]
    # Drop rows where all samples are 0 (TSSr operates on TSSrawMatrix which
    # only contains positive cells by construction)
    data = data.loc[~(data == 0).all(axis=1)]

    corr = data.corr(method="pearson")
    print("Correlation matrix (Pearson r):")
    print(corr.round(3).to_string())
    if output:
        corr.to_csv(output)
        typer.echo(f"Correlation matrix saved to {output}")
    if not plot:
        return

    _draw_pairs_plot(data, sample_cols, corr, plot_file, source.strip().lower())


def _draw_pairs_plot(data: pd.DataFrame, sample_cols, corr: pd.DataFrame,
                     plot_file: str, source: str):
    """R-pairs-style scatter + correlation panel (matches TSSr layout)."""
    n = len(sample_cols)
    fig, axes = plt.subplots(n, n, figsize=(2.5 * n, 2.5 * n), squeeze=False)
    for ax_row in axes:
        for ax in ax_row:
            ax.set_aspect("equal", adjustable="box")
            ax.set_xticks([]); ax.set_yticks([])

    use_log = (source == "raw")
    if not use_log:
        scatter_data = np.log10(data + 1)
        ax_lim = (0, float(scatter_data.values.max()) * 1.05)

    for i in range(n):
        for j in range(n):
            ax = axes[i, j]
            if i == j:
                # Diagonal — variable name (R pairs() default)
                ax.annotate(sample_cols[i], (0.5, 0.5), xycoords="axes fraction",
                            ha="center", va="center", fontsize=14, fontweight="bold")
                ax.set_frame_on(False)
            elif i > j:
                # Lower triangle — scatter
                if use_log:
                    x = data.iloc[:, j].astype(float)
                    y = data.iloc[:, i].astype(float)
                    mask = (x > 0) & (y > 0)
                    ax.scatter(x[mask], y[mask], s=1, color=TSSR_POINT_COLOR, alpha=0.6)
                    ax.set_xscale("log"); ax.set_yscale("log")
                else:
                    ax.scatter(scatter_data.iloc[:, j], scatter_data.iloc[:, i],
                               s=1, color=TSSR_POINT_COLOR, alpha=0.6)
                    ax.set_xlim(ax_lim); ax.set_ylim(ax_lim)
                if j == 0:
                    ax.set_ylabel(sample_cols[i])
                if i == n - 1:
                    ax.set_xlabel(sample_cols[j])
            else:
                # Upper triangle — Pearson r as text, font size scales with |r|
                r = corr.iloc[i, j]
                # TSSr: cex.cor = 0.8 / strwidth(txt); text size = cex.cor * r.
                # Approximation: scale font from ~8 (|r|≈0) to ~32 (|r|≈1)
                font_size = max(8.0, min(32.0, 8.0 + 24.0 * abs(r)))
                ax.annotate(f"{r:.2f}", (0.5, 0.5), xycoords="axes fraction",
                            ha="center", va="center",
                            fontsize=font_size, fontweight="bold")
                ax.set_frame_on(False)
    plt.tight_layout()
    plt.savefig(plot_file, dpi=150)
    typer.echo(f"Pairs plot saved to {plot_file}")
