#!/usr/bin/env python3
"""
Merge / normalize / filter — strict TSSr 0.99.6 parity.

Mirrors:
  TSSr::mergeSamples()   (MergingMethods.R)
  TSSr::normalizeTSS()   (NormalizationMethods.R)
  TSSr::filterTSS()      (FilteringMethods.R) — both 'poisson' and 'TPM' methods
"""

import typer
import pandas as pd
import numpy as np
from typing import Optional, List, Dict
from pathlib import Path
import logging

logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

__version__ = "0.14.0"

app = typer.Typer(help=f"Merge / normalize / filter TSS data (v{__version__})")

KEY_COLS = ['chr', 'pos', 'strand']


# -----------------------------------------------------------------------------
# Core algorithms (TSSr-parity)
# -----------------------------------------------------------------------------

def merge_samples_by_groups(df: pd.DataFrame,
                            sample_cols: List[str],
                            group_names: List[str],
                            merge_index: List[int]) -> pd.DataFrame:
    """TSSr-style mergeSamples: per-group rowSums of raw counts."""
    if len(merge_index) != len(sample_cols):
        raise ValueError(f"merge_index length ({len(merge_index)}) != sample count ({len(sample_cols)})")
    if len(set(merge_index)) != len(group_names):
        raise ValueError(f"unique mergeIndex count ({len(set(merge_index))}) != group_names count ({len(group_names)})")

    result = df[KEY_COLS].copy()
    for group_idx, group_name in enumerate(group_names, start=1):
        idxs = [i for i, m in enumerate(merge_index) if m == group_idx]
        cols = [sample_cols[i] for i in idxs]
        result[group_name] = df[cols].sum(axis=1)
    return result


def _is_normalized(df: pd.DataFrame, sample_cols: List[str]) -> bool:
    """TSSr's heuristic: data is normalized if any value is strictly in (0, 1)."""
    sub = df[sample_cols]
    return bool(((sub > 0) & (sub < 1)).any().any())


def normalize_to_tpm(df: pd.DataFrame,
                     sample_cols: List[str],
                     library_sizes: Optional[Dict[str, float]] = None) -> pd.DataFrame:
    """
    TSSr-style TPM normalization:  round(count / (lib_size / 1e6), 6)
    Sorts output by (strand, chr, pos) like TSSr's setorder().
    """
    if library_sizes is None:
        library_sizes = {col: float(df[col].sum()) for col in sample_cols}

    result = df.copy()
    for col in sample_cols:
        lib = library_sizes[col]
        if lib > 0:
            result[col] = (df[col].astype(float) / (lib / 1e6)).round(6)
        else:
            result[col] = 0.0
            logger.warning(f"Sample {col} has zero library size")
    return result.sort_values(["strand", "chr", "pos"]).reset_index(drop=True)


def filter_tss_by_poisson(df: pd.DataFrame,
                          sample_cols: List[str],
                          library_sizes: Dict[str, float],
                          genome_size: int,
                          p_val: float = 0.01,
                          normalization: bool = True) -> pd.DataFrame:
    """
    TSSr-style Poisson noise filter (FilteringFunctions.R::.filterWithPoisson):

        lambda = library_size / (genome_size * 2)
        cutoff = qpois(p_val, lambda, lower.tail=FALSE)
        cells with count < cutoff -> 0
        optionally normalize remaining to TPM (round 6)
        drop rows where all sample counts are 0
        sort by (strand, chr, pos)

    Requires raw integer counts. Raises if input looks already normalized.
    """
    from scipy.stats import poisson

    if _is_normalized(df, sample_cols):
        raise ValueError(
            "Poisson filter requires raw integer counts; input looks normalized. "
            "Run filter BEFORE normalize, or pass --no-normalize on prior step."
        )

    result = df.copy()
    for col in sample_cols:
        lib = library_sizes[col]
        lam = lib / (genome_size * 2)
        # R's qpois(p, lambda, lower.tail=FALSE) <=> scipy.stats.poisson.isf(p, lambda)
        cutoff = int(poisson.isf(p_val, lam))
        logger.info(f"  [{col}] lib={lib:,.0f}  lambda={lam:.6f}  cutoff={cutoff}")
        result.loc[result[col] < cutoff, col] = 0
        if normalization:
            result[col] = (result[col].astype(float) / (lib / 1e6)).round(6)

    # Drop all-zero rows
    keep = (result[sample_cols] > 0).any(axis=1)
    result = result[keep].copy()
    return result.sort_values(["strand", "chr", "pos"]).reset_index(drop=True)


def filter_tss_by_tpm(df: pd.DataFrame,
                      sample_cols: List[str],
                      tpm_low: float = 0.1) -> pd.DataFrame:
    """
    TSSr-style TPM filter:
        cells with TPM < tpm_low -> 0
        drop rows where all sample TPMs are 0
        sort by (strand, chr, pos)

    Requires normalized input.
    """
    if not _is_normalized(df, sample_cols):
        raise ValueError(
            "TPM filter requires normalized (TPM) data; input looks like raw counts. "
            "Run normalize before this filter."
        )

    result = df.copy()
    for col in sample_cols:
        result.loc[result[col] < tpm_low, col] = 0.0

    keep = (result[sample_cols] > 0).any(axis=1)
    result = result[keep].copy()
    return result.sort_values(["strand", "chr", "pos"]).reset_index(drop=True)


def compute_genome_size(fasta_path: str) -> int:
    """Sum of all chromosome lengths from FASTA (uses .fai index)."""
    import pysam
    fa = pysam.FastaFile(fasta_path)
    try:
        return sum(fa.get_reference_length(n) for n in fa.references)
    finally:
        fa.close()


def _sample_cols(df: pd.DataFrame) -> List[str]:
    return [c for c in df.columns if c not in KEY_COLS]


# -----------------------------------------------------------------------------
# CLI
# -----------------------------------------------------------------------------

@app.command("merge")
def merge_command(
    input_file: Path = typer.Option(..., "-i", "--input", help="Raw TSS table from tssCalling"),
    output_file: Path = typer.Option(..., "-o", "--output", help="Output merged TSS table"),
    group_names: str = typer.Option(..., "-g", "--groups", help='Group names (e.g. "control treat")'),
    merge_index: str = typer.Option(..., "-m", "--merge-index", help='1-indexed per-sample group (e.g. "1 1 2 2")'),
):
    """Merge biological replicates by group (TSSr mergeSamples)."""
    df = pd.read_csv(input_file, sep="\t")
    sample_cols = _sample_cols(df)
    indices = [int(x) for x in merge_index.split()]
    groups = group_names.split()
    logger.info(f"Samples: {sample_cols}  Groups: {groups}  mergeIndex: {indices}")
    out = merge_samples_by_groups(df, sample_cols, groups, indices)
    out.to_csv(output_file, sep="\t", index=False)
    for col in _sample_cols(out):
        logger.info(f"  library size [{col}]: {int(out[col].sum()):,}")
    logger.info(f"Saved {len(out):,} rows to {output_file}")


@app.command("normalize")
def normalize_command(
    input_file: Path = typer.Option(..., "-i", "--input", help="Merged (raw count) TSS table"),
    output_file: Path = typer.Option(..., "-o", "--output", help="Output normalized TSS table"),
):
    """TPM normalize (TSSr normalizeTSS): round(count/(lib_size/1e6), 6)."""
    df = pd.read_csv(input_file, sep="\t")
    sample_cols = _sample_cols(df)
    out = normalize_to_tpm(df, sample_cols)
    out.to_csv(output_file, sep="\t", index=False)
    logger.info(f"Saved {len(out):,} rows to {output_file}")


@app.command("filter")
def filter_command(
    input_file: Path = typer.Option(..., "-i", "--input", help="TSS table to filter"),
    output_file: Path = typer.Option(..., "-o", "--output", help="Output filtered TSS table"),
    method: str = typer.Option("poisson", "--method", help="'poisson' (default) or 'TPM'"),
    p_val: float = typer.Option(0.01, "--p-val", help="Poisson p-value threshold (default 0.01)"),
    tpm_low: float = typer.Option(0.1, "--tpm-low", help="TPM threshold (default 0.1)"),
    normalization: bool = typer.Option(True, "--normalization/--no-normalization",
                                       help="For poisson: normalize to TPM after filtering"),
    reference: Optional[Path] = typer.Option(None, "-r", "--reference",
                                             help="Reference FASTA (needed for poisson — for genome size)"),
    genome_size: Optional[int] = typer.Option(None, "--genome-size",
                                              help="Genome size in bp (alternative to --reference)"),
):
    """
    Filter TSS data (TSSr filterTSS).

    'poisson' requires raw counts and either --reference or --genome-size.
    'TPM' requires normalized input.
    """
    df = pd.read_csv(input_file, sep="\t")
    sample_cols = _sample_cols(df)

    method_lc = method.strip().lower()
    if method_lc == "poisson":
        if genome_size is None:
            if reference is None:
                raise typer.BadParameter("poisson filter needs --reference (FASTA) or --genome-size")
            genome_size = compute_genome_size(str(reference))
            logger.info(f"genome_size (from FASTA): {genome_size:,}")
        lib_sizes = {col: float(df[col].sum()) for col in sample_cols}
        out = filter_tss_by_poisson(df, sample_cols, lib_sizes, genome_size,
                                    p_val=p_val, normalization=normalization)
    elif method_lc == "tpm":
        out = filter_tss_by_tpm(df, sample_cols, tpm_low=tpm_low)
    else:
        raise typer.BadParameter(f"Unknown method: {method!r}. Use 'poisson' or 'TPM'.")

    out.to_csv(output_file, sep="\t", index=False)
    logger.info(f"Filtered {len(df):,} → {len(out):,} rows. Saved to {output_file}")


@app.command("process")
def process_command(
    input_file: Path = typer.Option(..., "-i", "--input", help="Raw TSS table from tssCalling"),
    output_file: Path = typer.Option(..., "-o", "--output", help="Output processed TSS table"),
    group_names: Optional[str] = typer.Option(None, "-g", "--groups"),
    merge_index: Optional[str] = typer.Option(None, "-m", "--merge-index"),
    filter_method: Optional[str] = typer.Option(
        None, "--filter",
        help="'poisson' or 'TPM' (omit to skip filtering)"),
    p_val: float = typer.Option(0.01, "--p-val"),
    tpm_low: float = typer.Option(0.1, "--tpm-low"),
    normalize: bool = typer.Option(True, "--normalize/--no-normalize",
                                   help="Normalize to TPM. With --filter poisson, applied AFTER filter."),
    reference: Optional[Path] = typer.Option(None, "-r", "--reference",
                                             help="Reference FASTA (for poisson filter genome size)"),
    genome_size: Optional[int] = typer.Option(None, "--genome-size"),
):
    """
    One-shot pipeline: merge → (filter poisson | normalize | filter TPM).

    Matches TSSr workflow:
      - tssCalling (upstream) → process: groups via -g/-m
      - For poisson filter: filter operates on raw counts, then normalizes.
      - For TPM filter: normalize first, then filter.
    """
    df = pd.read_csv(input_file, sep="\t")
    sample_cols = _sample_cols(df)
    logger.info(f"Input: {len(df):,} rows, samples: {sample_cols}")

    # Step 1: merge
    if group_names and merge_index:
        groups = group_names.split()
        indices = [int(x) for x in merge_index.split()]
        df = merge_samples_by_groups(df, sample_cols, groups, indices)
        sample_cols = _sample_cols(df)
        logger.info(f"After merge: {sample_cols}")

    lib_sizes = {col: float(df[col].sum()) for col in sample_cols}

    method_lc = filter_method.strip().lower() if filter_method else None

    if method_lc == "poisson":
        if genome_size is None:
            if reference is None:
                raise typer.BadParameter("poisson filter needs --reference or --genome-size")
            genome_size = compute_genome_size(str(reference))
            logger.info(f"genome_size: {genome_size:,}")
        df = filter_tss_by_poisson(df, sample_cols, lib_sizes, genome_size,
                                   p_val=p_val, normalization=normalize)
    elif method_lc == "tpm":
        # TSSr semantics: TPM filter requires normalized input → normalize first.
        if normalize:
            df = normalize_to_tpm(df, sample_cols, library_sizes=lib_sizes)
        df = filter_tss_by_tpm(df, sample_cols, tpm_low=tpm_low)
    elif method_lc is None:
        if normalize:
            df = normalize_to_tpm(df, sample_cols, library_sizes=lib_sizes)
    else:
        raise typer.BadParameter(f"Unknown filter method: {filter_method!r}")

    df.to_csv(output_file, sep="\t", index=False)
    logger.info(f"Final: {len(df):,} rows saved to {output_file}")


if __name__ == "__main__":
    app()
