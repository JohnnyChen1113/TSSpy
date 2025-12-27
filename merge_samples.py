#!/usr/bin/env python3
"""
Merge Samples - Merge biological replicates and normalize TSS data
"""

import typer
import pandas as pd
import numpy as np
from typing import Optional, List
from pathlib import Path
import logging

logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

__version__ = "0.1.0"

app = typer.Typer(help=f"Merge samples and normalize TSS data (v{__version__})")


def merge_samples_by_groups(df: pd.DataFrame,
                            sample_cols: List[str],
                            group_names: List[str],
                            merge_index: List[int]) -> pd.DataFrame:
    """
    Merge samples according to group assignments.

    Args:
        df: TSS table with chr, pos, strand, and sample columns
        sample_cols: List of sample column names
        group_names: Names for merged groups
        merge_index: Group assignment for each sample (1-indexed)

    Returns:
        DataFrame with merged sample columns
    """
    result = df[['chr', 'pos', 'strand']].copy()

    # Create merged columns
    for group_idx, group_name in enumerate(group_names, start=1):
        # Find samples belonging to this group
        sample_indices = [i for i, idx in enumerate(merge_index) if idx == group_idx]
        if not sample_indices:
            logger.warning(f"No samples assigned to group {group_idx} ({group_name})")
            continue

        samples_to_merge = [sample_cols[i] for i in sample_indices]
        logger.info(f"Merging samples {samples_to_merge} into '{group_name}'")

        # Sum the counts from all samples in this group
        result[group_name] = df[samples_to_merge].sum(axis=1)

    return result


def normalize_to_tpm(df: pd.DataFrame, sample_cols: List[str]) -> pd.DataFrame:
    """
    Normalize raw counts to TPM (Tags Per Million).

    Args:
        df: TSS table
        sample_cols: List of sample columns to normalize

    Returns:
        DataFrame with TPM values
    """
    result = df.copy()

    for col in sample_cols:
        total = df[col].sum()
        if total > 0:
            result[col] = df[col] / total * 1e6
        else:
            result[col] = 0.0
            logger.warning(f"Sample {col} has zero total counts")

    return result


def filter_tss_by_tpm(df: pd.DataFrame,
                      sample_cols: List[str],
                      threshold: float = 0.1) -> pd.DataFrame:
    """
    Filter TSS positions by TPM threshold.

    A position is kept if at least one sample has TPM >= threshold.

    Args:
        df: TSS table with TPM values
        sample_cols: Sample columns to check
        threshold: Minimum TPM value

    Returns:
        Filtered DataFrame
    """
    mask = (df[sample_cols] >= threshold).any(axis=1)
    filtered = df[mask].copy()
    logger.info(f"Filtered from {len(df)} to {len(filtered)} positions (TPM >= {threshold})")
    return filtered


def filter_tss_by_poisson(df: pd.DataFrame,
                          sample_cols: List[str],
                          library_sizes: dict,
                          pvalue_threshold: float = 0.01) -> pd.DataFrame:
    """
    Filter TSS positions using Poisson distribution.

    Only keep TSS positions where the observed count is significantly
    higher than expected by chance.

    Args:
        df: TSS table with raw counts
        sample_cols: Sample columns to filter
        library_sizes: Dictionary mapping sample names to library sizes
        pvalue_threshold: P-value threshold

    Returns:
        Filtered DataFrame
    """
    from scipy.stats import poisson

    # Calculate expected count per position based on library size
    # Null hypothesis: reads are uniformly distributed across genome
    # For simplicity, we use a Poisson test comparing to the average

    keep_mask = pd.Series([False] * len(df), index=df.index)

    for col in sample_cols:
        lib_size = library_sizes.get(col, df[col].sum())
        # Expected count per position (simplified)
        n_positions = len(df)
        expected_lambda = lib_size / n_positions if n_positions > 0 else 1

        # Calculate p-value for each position
        for idx, count in df[col].items():
            if count > 0:
                # P(X >= count) for Poisson distribution
                pval = 1 - poisson.cdf(count - 1, expected_lambda)
                if pval < pvalue_threshold:
                    keep_mask[idx] = True

    filtered = df[keep_mask].copy()
    logger.info(f"Poisson filter: {len(df)} -> {len(filtered)} positions (p < {pvalue_threshold})")
    return filtered


@app.command("merge")
def merge_command(
    input_file: Path = typer.Option(
        ..., "-i", "--input",
        help="Input TSS table file"
    ),
    output_file: Path = typer.Option(
        ..., "-o", "--output",
        help="Output merged TSS table file"
    ),
    group_names: str = typer.Option(
        ..., "-g", "--groups",
        help="Group names (space-separated)"
    ),
    merge_index: str = typer.Option(
        ..., "-m", "--merge-index",
        help="Group assignment for each sample (space-separated integers, 1-indexed)"
    ),
):
    """
    Merge biological replicates into groups.

    Example:
        tsspy mergeSamples merge -i raw.tsv -o merged.tsv -g "control treat" -m "1 1 2 2"

        This merges samples 1,2 into "control" and samples 3,4 into "treat"
    """
    logger.info(f"Reading input file: {input_file}")
    df = pd.read_csv(input_file, sep='\t')

    # Parse arguments
    groups = group_names.split()
    indices = [int(x) for x in merge_index.split()]

    # Get sample columns (everything except chr, pos, strand)
    sample_cols = [c for c in df.columns if c not in ['chr', 'pos', 'strand']]

    if len(indices) != len(sample_cols):
        raise typer.BadParameter(
            f"Number of merge indices ({len(indices)}) must match number of samples ({len(sample_cols)})"
        )

    logger.info(f"Samples: {sample_cols}")
    logger.info(f"Groups: {groups}")
    logger.info(f"Merge index: {indices}")

    # Merge samples
    result = merge_samples_by_groups(df, sample_cols, groups, indices)

    # Save
    result.to_csv(output_file, sep='\t', index=False)
    logger.info(f"Saved merged table to {output_file}")

    # Report library sizes
    for col in [c for c in result.columns if c not in ['chr', 'pos', 'strand']]:
        lib_size = result[col].sum()
        logger.info(f"Library size for {col}: {lib_size:,}")


@app.command("normalize")
def normalize_command(
    input_file: Path = typer.Option(
        ..., "-i", "--input",
        help="Input TSS table file"
    ),
    output_file: Path = typer.Option(
        ..., "-o", "--output",
        help="Output normalized TSS table file"
    ),
):
    """
    Normalize TSS counts to TPM (Tags Per Million).

    Example:
        tsspy mergeSamples normalize -i merged.tsv -o normalized.tsv
    """
    logger.info(f"Reading input file: {input_file}")
    df = pd.read_csv(input_file, sep='\t')

    sample_cols = [c for c in df.columns if c not in ['chr', 'pos', 'strand']]

    # Report original library sizes
    for col in sample_cols:
        lib_size = df[col].sum()
        logger.info(f"Library size for {col}: {lib_size:,}")

    # Normalize
    result = normalize_to_tpm(df, sample_cols)

    # Save
    result.to_csv(output_file, sep='\t', index=False)
    logger.info(f"Saved normalized table to {output_file}")


@app.command("filter")
def filter_command(
    input_file: Path = typer.Option(
        ..., "-i", "--input",
        help="Input TSS table file"
    ),
    output_file: Path = typer.Option(
        ..., "-o", "--output",
        help="Output filtered TSS table file"
    ),
    method: str = typer.Option(
        "tpm", "--method",
        help="Filter method: 'tpm' or 'poisson'"
    ),
    threshold: float = typer.Option(
        0.1, "--threshold",
        help="TPM threshold (for tpm method) or p-value threshold (for poisson)"
    ),
):
    """
    Filter TSS positions by expression threshold.

    Example:
        tsspy mergeSamples filter -i normalized.tsv -o filtered.tsv --method tpm --threshold 0.1
    """
    logger.info(f"Reading input file: {input_file}")
    df = pd.read_csv(input_file, sep='\t')

    sample_cols = [c for c in df.columns if c not in ['chr', 'pos', 'strand']]

    if method.lower() == 'tpm':
        result = filter_tss_by_tpm(df, sample_cols, threshold)
    elif method.lower() == 'poisson':
        library_sizes = {col: df[col].sum() for col in sample_cols}
        result = filter_tss_by_poisson(df, sample_cols, library_sizes, threshold)
    else:
        raise typer.BadParameter(f"Unknown filter method: {method}")

    # Save
    result.to_csv(output_file, sep='\t', index=False)
    logger.info(f"Saved filtered table to {output_file}")


@app.command("process")
def process_command(
    input_file: Path = typer.Option(
        ..., "-i", "--input",
        help="Input raw TSS table file"
    ),
    output_file: Path = typer.Option(
        ..., "-o", "--output",
        help="Output processed TSS table file"
    ),
    group_names: Optional[str] = typer.Option(
        None, "-g", "--groups",
        help="Group names (space-separated). If not provided, no merging is done."
    ),
    merge_index: Optional[str] = typer.Option(
        None, "-m", "--merge-index",
        help="Group assignment for each sample (space-separated integers)"
    ),
    normalize: bool = typer.Option(
        True, "--normalize/--no-normalize",
        help="Normalize to TPM"
    ),
    filter_method: Optional[str] = typer.Option(
        None, "--filter",
        help="Filter method: 'tpm' or 'poisson'"
    ),
    filter_threshold: float = typer.Option(
        0.1, "--filter-threshold",
        help="Filter threshold value"
    ),
):
    """
    One-step processing: merge, normalize, and filter TSS data.

    Example:
        tsspy mergeSamples process -i raw.tsv -o processed.tsv \\
            -g "control treat" -m "1 1 2 2" \\
            --normalize --filter tpm --filter-threshold 0.1
    """
    logger.info(f"Reading input file: {input_file}")
    df = pd.read_csv(input_file, sep='\t')

    sample_cols = [c for c in df.columns if c not in ['chr', 'pos', 'strand']]
    logger.info(f"Input samples: {sample_cols}")

    # Step 1: Merge (if specified)
    if group_names and merge_index:
        groups = group_names.split()
        indices = [int(x) for x in merge_index.split()]
        if len(indices) != len(sample_cols):
            raise typer.BadParameter(
                f"Number of merge indices ({len(indices)}) must match number of samples ({len(sample_cols)})"
            )
        df = merge_samples_by_groups(df, sample_cols, groups, indices)
        sample_cols = [c for c in df.columns if c not in ['chr', 'pos', 'strand']]
        logger.info(f"After merge: {sample_cols}")

    # Step 2: Normalize (if specified)
    if normalize:
        df = normalize_to_tpm(df, sample_cols)
        logger.info("Normalized to TPM")

    # Step 3: Filter (if specified)
    if filter_method:
        if filter_method.lower() == 'tpm':
            df = filter_tss_by_tpm(df, sample_cols, filter_threshold)
        elif filter_method.lower() == 'poisson':
            library_sizes = {col: df[col].sum() for col in sample_cols}
            df = filter_tss_by_poisson(df, sample_cols, library_sizes, filter_threshold)

    # Save
    df.to_csv(output_file, sep='\t', index=False)
    logger.info(f"Saved processed table to {output_file}")
    logger.info(f"Final: {len(df)} TSS positions, {len(sample_cols)} samples")


if __name__ == '__main__':
    app()
