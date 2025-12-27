#!/usr/bin/env python3
"""
TSSpy: Python CLI for TSS Analysis

A comprehensive command-line tool for transcription start site (TSS) data analysis.
Inspired by TSSr (R/Bioconductor) but implemented in pure Python without BSgenome dependency.

Main features:
- TSS calling from BAM files with reference-based G mismatch removal
- Sample merging and normalization (TPM)
- TSS clustering to infer core promoters
- Consensus clustering across samples
- Promoter shape analysis (PSS/SI scores)
- Gene assignment from GTF/GFF annotations
- BigWig/BedGraph export for visualization
- Correlation analysis and plotting

Usage:
    tsspy <command> [options]

Commands:
    tssCalling       - Extract TSS from BAM files
    mergeSamples     - Merge samples and normalize data
    clustering       - Cluster TSS to infer core promoters
    consensusCluster - Create consensus clusters across samples
    shapeCluster     - Calculate promoter shape scores (PSS/SI)
    geneAssign       - Assign clusters to genes
    bigwig           - Generate BigWig/BedGraph files
    correlation      - Calculate sample correlations
    plot             - Generate visualization plots
"""

import typer

from TSSpy import tss_calling
from TSSpy import clustering
from TSSpy import gene_assign
from TSSpy import plot
from TSSpy import bigwig
from TSSpy import merge_samples
from TSSpy import consensus_cluster
from TSSpy import shape_cluster
from TSSpy.correlation import correlation

__version__ = "0.3.0"

# Create main app
app = typer.Typer(
    name="tsspy",
    help=f"TSSpy: Python CLI for TSS analysis (v{__version__})",
    add_completion=False,
)

# Register all subcommands
app.add_typer(tss_calling.app, name="tssCalling", help="Extract TSS from BAM files")
app.add_typer(merge_samples.app, name="mergeSamples", help="Merge samples and normalize")
app.add_typer(clustering.app, name="clustering", help="Cluster TSS to infer core promoters")
app.add_typer(consensus_cluster.app, name="consensusCluster", help="Consensus clustering across samples")
app.add_typer(shape_cluster.app, name="shapeCluster", help="Calculate promoter shape scores")
app.add_typer(gene_assign.app, name="geneAssign", help="Assign clusters to genes")
app.add_typer(bigwig.app, name="bigwig", help="Generate BigWig/BedGraph files")
app.add_typer(plot.app, name="plot", help="Generate visualization plots")

# Register correlation as a direct command
app.command(name="correlation", help="Calculate sample correlations")(correlation)


@app.command()
def version():
    """Show version information."""
    typer.echo(f"TSSpy version {__version__}")
    typer.echo("A Python CLI for TSS analysis")
    typer.echo("https://github.com/JohnnyChen1113/TSSpy")


@app.callback()
def main_callback(
    ctx: typer.Context,
    verbose: bool = typer.Option(False, "--verbose", "-v", help="Enable verbose output"),
):
    """
    TSSpy: Python CLI for TSS Analysis

    A comprehensive tool for analyzing transcription start site (TSS) data.
    Inspired by TSSr but implemented in pure Python.
    """
    if verbose:
        import logging
        logging.basicConfig(level=logging.DEBUG)


def main():
    """Entry point for the CLI."""
    app()


if __name__ == "__main__":
    main()
