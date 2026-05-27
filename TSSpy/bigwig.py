#!/usr/bin/env python3
"""
Export TSS table to bedGraph / BigWig — strict parity with TSSr 0.99.6
exportTSStoBedgraph(format = "bedGraph" | "BigWig").

Per-sample × per-strand output files:
   <sample>.TSS.<data>.plus.bedGraph
   <sample>.TSS.<data>.minus.bedGraph
   <sample>.TSS.<data>.plus.BigWig
   <sample>.TSS.<data>.minus.BigWig

bedGraph format: tab-separated  chr  start_0based  end  value
   (no track header — rtracklayer's default export doesn't write one)
Minus-strand scores are negated (TSSr convention).

Input: TSS table (chr, pos, strand, sample columns). Pass the merged/
filtered table from mergeSamples/filterTSS to reproduce TSSr's
data="processed" path. No additional normalization is done here.
"""

from __future__ import annotations
import typer
import pandas as pd
import numpy as np
import pyBigWig
import pysam
from pathlib import Path
from typing import List, Optional, Dict
import logging

logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

__version__ = "0.19.0"
app = typer.Typer(help=f"Export TSS table to bedGraph / BigWig (v{__version__})")

KEY_COLS = ['chr', 'pos', 'strand']


def _format_value(v: float) -> str:
    """Match rtracklayer's default float formatting (Python repr / shortest roundtrip)."""
    # Integer-valued floats are written without decimal in TSSr (e.g. "2", not "2.0")
    if v == int(v):
        return str(int(v))
    return repr(v)


def _read_chrom_sizes(fasta_path: Optional[str], bam_path: Optional[str]) -> Dict[str, int]:
    """Build chrom->size from FASTA .fai or from a BAM header."""
    if fasta_path is not None:
        fa = pysam.FastaFile(fasta_path)
        try:
            return {n: fa.get_reference_length(n) for n in fa.references}
        finally:
            fa.close()
    if bam_path is not None:
        with pysam.AlignmentFile(bam_path, "rb") as bam:
            return dict(zip(bam.references, bam.lengths))
    raise ValueError("Need --reference FASTA or --bam for chromosome sizes.")


def write_bedgraph(df: pd.DataFrame, sample: str, strand: str, out_path: str) -> int:
    """Write TSSr-style bedGraph for one sample × strand. Returns row count."""
    sub = df[(df['strand'] == strand)][['chr', 'pos', sample]].copy()
    sub = sub[sub[sample] != 0]
    if sub.empty:
        Path(out_path).touch()
        return 0
    sub = sub.sort_values(['chr', 'pos'], kind='stable').reset_index(drop=True)
    with open(out_path, 'w') as f:
        for chrom, pos, val in zip(sub['chr'].values, sub['pos'].values, sub[sample].values):
            f.write(f"{chrom}\t{int(pos) - 1}\t{int(pos)}\t{_format_value(float(val))}\n")
    return len(sub)


def write_bigwig(df: pd.DataFrame, sample: str, strand: str,
                 chrom_sizes: Dict[str, int], out_path: str) -> int:
    """Write TSSr-style BigWig for one sample × strand. Returns row count.

    Only chromosomes that have data appear in the BigWig header — matches
    TSSr's `seqlengths(temp.p) <- seqlengths(Genome)[... %in% temp.p@seqnames]`.
    """
    sub = df[(df['strand'] == strand)][['chr', 'pos', sample]].copy()
    sub = sub[sub[sample] != 0]
    # Restrict to chromosomes in the size table; preserve their order
    sub = sub[sub['chr'].isin(chrom_sizes)]
    # Header only includes chromosomes that actually have data, in FASTA order
    present = [c for c in chrom_sizes.keys() if c in set(sub['chr'].unique())]
    chr_order = {c: i for i, c in enumerate(present)}
    sub['_o'] = sub['chr'].map(chr_order)
    sub = sub.sort_values(['_o', 'pos'], kind='stable').reset_index(drop=True)

    bw = pyBigWig.open(out_path, "w")
    bw.addHeader([(c, chrom_sizes[c]) for c in present])
    if not sub.empty:
        for chrom, grp in sub.groupby('chr', sort=False):
            chroms = [chrom] * len(grp)
            starts = (grp['pos'].astype(int) - 1).tolist()
            ends = grp['pos'].astype(int).tolist()
            values = grp[sample].astype(float).tolist()
            bw.addEntries(chroms, starts, ends=ends, values=values)
    bw.close()
    return len(sub)


@app.callback(invoke_without_command=True)
def main(
    ctx: typer.Context,
    input_file: Optional[Path] = typer.Option(None, "-i", "--input",
                                              help="TSS table (chr, pos, strand, <sample>...)"),
    output_prefix: Optional[str] = typer.Option(None, "-o", "--output-prefix",
                                                help="Output prefix (e.g. 'YPD.TSS.processed' → '<prefix>.plus.bedGraph' etc.) — "
                                                     "if multiple samples, use --batch to produce one prefix per sample"),
    samples: Optional[str] = typer.Option(None, "-s", "--samples",
                                          help="Comma-separated sample columns. Default: all non-key columns."),
    fmt: str = typer.Option("bedGraph", "--format",
                            help="'bedGraph' or 'BigWig'"),
    reference: Optional[Path] = typer.Option(None, "-r", "--reference",
                                             help="Reference FASTA (for BigWig chrom sizes)"),
    bam_for_sizes: Optional[Path] = typer.Option(None, "--bam",
                                                 help="BAM file (alternative source of chrom sizes)"),
    data_label: str = typer.Option("processed", "--data-label",
                                   help="Filename infix matching TSSr's `data` arg ('processed' or 'raw')"),
    batch: bool = typer.Option(False, "--batch",
                               help="Per-sample filenames matching TSSr: <sample>.TSS.<data>.{plus,minus}.{bedGraph,BigWig}"),
    version: bool = typer.Option(False, "--version"),
):
    """Export TSS table to bedGraph / BigWig (TSSr exportTSStoBedgraph parity)."""
    if version:
        typer.echo(f"bigwig.py v{__version__}")
        raise typer.Exit(0)
    if input_file is None:
        typer.echo(ctx.get_help())
        raise typer.Exit(0)

    df = pd.read_csv(input_file, sep="\t")
    sample_cols = ([s.strip() for s in samples.split(',')] if samples
                   else [c for c in df.columns if c not in KEY_COLS])
    fmt_lc = fmt.strip().lower()
    if fmt_lc not in {"bedgraph", "bigwig"}:
        raise typer.BadParameter("--format must be 'bedGraph' or 'BigWig'")

    chrom_sizes: Optional[Dict[str, int]] = None
    if fmt_lc == "bigwig":
        chrom_sizes = _read_chrom_sizes(str(reference) if reference else None,
                                         str(bam_for_sizes) if bam_for_sizes else None)

    # Negate minus-strand values once at the dataframe level (matches TSSr's
    # `temp[, score := score * -1]` for minus).
    for sc in sample_cols:
        # Make a per-sample copy so per-sample writes don't double-negate.
        pass  # handled below per write

    for sc in sample_cols:
        ext = "bedGraph" if fmt_lc == "bedgraph" else "BigWig"
        if batch:
            base = f"{sc}.TSS.{data_label}"
        else:
            base = (output_prefix or sc)
        plus_path = f"{base}.plus.{ext}"
        minus_path = f"{base}.minus.{ext}"

        # For minus strand: write negated values
        df_for_sample = df[KEY_COLS + [sc]].copy()
        df_for_sample.loc[df_for_sample['strand'] == '-', sc] = (
            -df_for_sample.loc[df_for_sample['strand'] == '-', sc].astype(float))

        if fmt_lc == "bedgraph":
            n_p = write_bedgraph(df_for_sample, sc, '+', plus_path)
            n_m = write_bedgraph(df_for_sample, sc, '-', minus_path)
        else:
            n_p = write_bigwig(df_for_sample, sc, '+', chrom_sizes, plus_path)
            n_m = write_bigwig(df_for_sample, sc, '-', chrom_sizes, minus_path)
        logger.info(f"{sc}: +={n_p}  -={n_m}  -> {plus_path}, {minus_path}")


if __name__ == '__main__':
    app()
