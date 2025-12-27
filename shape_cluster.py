#!/usr/bin/env python3
"""
Shape Cluster - Calculate core promoter shape scores (PSS and SI)

PSS (Promoter Shape Score): Lu and Lin 2019
SI (Shape Index): Hoskins et al. 2011
"""

import typer
import pandas as pd
import numpy as np
from typing import Optional, List
from pathlib import Path
import logging
from multiprocessing import Pool, cpu_count

logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

__version__ = "0.1.0"

app = typer.Typer(help=f"Calculate core promoter shape scores (v{__version__})")


def calculate_pss(tags: np.ndarray, interquantile_width: int) -> float:
    """
    Calculate Promoter Shape Score (PSS).

    PSS = -sum(p_i * log2(p_i)) * log2(IQW)

    Where:
        - p_i = tags_i / total_tags (probability of each TSS)
        - IQW = interquantile width (q_0.9 - q_0.1 + 1)

    PSS integrates both the distribution entropy and the width of the promoter.
    Lower PSS = sharper promoter (more focused TSS signal)
    PSS = 0 for singletons (single TSS position)

    Args:
        tags: Array of tag counts for TSS positions within the cluster
        interquantile_width: Width of the interquantile region (q_0.9 - q_0.1 + 1)

    Returns:
        Promoter Shape Score
    """
    if len(tags) == 0 or interquantile_width <= 0:
        return 0.0

    total = tags.sum()
    if total == 0:
        return 0.0

    # Handle singleton case
    if len(tags) == 1 or interquantile_width == 1:
        return 0.0

    # Calculate probabilities
    probs = tags / total

    # Calculate entropy: -sum(p * log2(p))
    # Avoid log(0) by filtering out zero probabilities
    nonzero_mask = probs > 0
    if not nonzero_mask.any():
        return 0.0

    entropy = -np.sum(probs[nonzero_mask] * np.log2(probs[nonzero_mask]))

    # PSS = entropy * log2(IQW)
    pss = entropy * np.log2(interquantile_width)

    return pss


def calculate_si(tags: np.ndarray) -> float:
    """
    Calculate Shape Index (SI).

    SI = 2 + sum(p_i * log2(p_i))

    Where p_i = tags_i / total_tags

    SI is based on Shannon entropy.
    Higher SI = sharper promoter
    SI = 2 for singletons (maximum sharpness)
    SI approaches 0 for completely uniform distribution

    Args:
        tags: Array of tag counts for TSS positions within the cluster

    Returns:
        Shape Index
    """
    if len(tags) == 0:
        return 0.0

    total = tags.sum()
    if total == 0:
        return 0.0

    # Singleton case
    if len(tags) == 1:
        return 2.0

    # Calculate probabilities
    probs = tags / total

    # Calculate entropy term: sum(p * log2(p))
    # Note: this is negative of standard entropy
    nonzero_mask = probs > 0
    if not nonzero_mask.any():
        return 2.0

    entropy_term = np.sum(probs[nonzero_mask] * np.log2(probs[nonzero_mask]))

    # SI = 2 + entropy_term (entropy_term is negative, so this reduces SI)
    si = 2.0 + entropy_term

    return si


def get_tss_in_cluster(tss_df: pd.DataFrame, cluster: pd.Series,
                       sample_col: str, use_interquantile: bool = True) -> np.ndarray:
    """
    Extract TSS tag counts within a cluster's region.

    Args:
        tss_df: TSS table with chr, pos, strand, and sample columns
        cluster: Series with cluster info (chr, strand, start, end, q_0.1, q_0.9)
        sample_col: Name of sample column to use
        use_interquantile: If True, use q_0.1 to q_0.9 region; otherwise use start to end

    Returns:
        Array of tag counts for TSS positions in the region
    """
    chr_name = cluster['chr']
    strand = cluster['strand']

    if use_interquantile and 'q_0.1' in cluster.index and 'q_0.9' in cluster.index:
        start = cluster['q_0.1']
        end = cluster['q_0.9']
    else:
        start = cluster['start']
        end = cluster['end']

    # Filter TSS positions within the region
    mask = (
        (tss_df['chr'] == chr_name) &
        (tss_df['strand'] == strand) &
        (tss_df['pos'] >= start) &
        (tss_df['pos'] <= end)
    )

    return tss_df.loc[mask, sample_col].values


def calculate_shape_for_cluster(args):
    """
    Calculate shape score for a single cluster.

    Args:
        args: Tuple of (cluster_dict, tss_df, sample_col, method)

    Returns:
        Dictionary with cluster info and shape score
    """
    cluster_dict, tss_df, sample_col, method = args

    cluster = pd.Series(cluster_dict)

    # Get TSS tags in the interquantile region
    tags = get_tss_in_cluster(tss_df, cluster, sample_col, use_interquantile=True)

    # Calculate interquantile width
    if 'q_0.1' in cluster.index and 'q_0.9' in cluster.index:
        iqw = int(cluster['q_0.9'] - cluster['q_0.1'] + 1)
    elif 'interquantile_width' in cluster.index:
        iqw = int(cluster['interquantile_width'])
    else:
        iqw = int(cluster['end'] - cluster['start'] + 1)

    # Calculate shape score
    if method.upper() == 'PSS':
        shape_score = calculate_pss(tags, iqw)
    elif method.upper() == 'SI':
        shape_score = calculate_si(tags)
    else:
        shape_score = np.nan

    result = cluster_dict.copy()
    result['shape_score'] = shape_score

    return result


def calculate_shape_scores(cluster_df: pd.DataFrame,
                            tss_df: pd.DataFrame,
                            sample_col: str,
                            method: str = 'PSS',
                            n_processes: int = 1) -> pd.DataFrame:
    """
    Calculate shape scores for all clusters.

    Args:
        cluster_df: DataFrame with cluster information
        tss_df: TSS table
        sample_col: Sample column to use for tag counts
        method: 'PSS' or 'SI'
        n_processes: Number of processes for parallel computation

    Returns:
        DataFrame with shape scores added
    """
    # Prepare arguments for each cluster
    cluster_dicts = cluster_df.to_dict('records')
    args_list = [(c, tss_df, sample_col, method) for c in cluster_dicts]

    # Calculate shape scores
    if n_processes > 1 and len(args_list) > 100:
        with Pool(processes=n_processes) as pool:
            results = pool.map(calculate_shape_for_cluster, args_list)
    else:
        results = [calculate_shape_for_cluster(args) for args in args_list]

    return pd.DataFrame(results)


@app.command("calculate")
def calculate_command(
    cluster_file: Path = typer.Option(
        ..., "-c", "--clusters",
        help="Input cluster file (from clustering or consensusCluster)"
    ),
    tss_file: Path = typer.Option(
        ..., "-t", "--tss",
        help="Input TSS table (processed/normalized)"
    ),
    output_file: Path = typer.Option(
        ..., "-o", "--output",
        help="Output file with shape scores"
    ),
    sample: str = typer.Option(
        ..., "-s", "--sample",
        help="Sample column to use for calculation"
    ),
    method: str = typer.Option(
        "PSS", "-m", "--method",
        help="Shape calculation method: PSS or SI"
    ),
    processes: int = typer.Option(
        1, "-p", "--processes",
        help="Number of processes"
    ),
):
    """
    Calculate promoter shape scores for clusters.

    PSS (Promoter Shape Score):
        - Combines entropy and interquantile width
        - Lower value = sharper promoter
        - PSS = 0 for singletons

    SI (Shape Index):
        - Based on Shannon entropy
        - Higher value = sharper promoter
        - SI = 2 for singletons

    Example:
        tsspy shapeCluster calculate -c clusters.tsv -t tss.tsv -o shape.tsv -s control -m PSS
    """
    logger.info(f"Loading cluster file: {cluster_file}")
    cluster_df = pd.read_csv(cluster_file, sep='\t')

    logger.info(f"Loading TSS file: {tss_file}")
    tss_df = pd.read_csv(tss_file, sep='\t')

    # Check if sample column exists
    if sample not in tss_df.columns:
        available = [c for c in tss_df.columns if c not in ['chr', 'pos', 'strand']]
        raise typer.BadParameter(f"Sample '{sample}' not found. Available: {available}")

    logger.info(f"Calculating {method} shape scores for sample '{sample}'")
    logger.info(f"Number of clusters: {len(cluster_df)}")

    result_df = calculate_shape_scores(
        cluster_df, tss_df, sample, method=method, n_processes=processes
    )

    # Save results
    result_df.to_csv(output_file, sep='\t', index=False)
    logger.info(f"Saved results to {output_file}")

    # Report statistics
    scores = result_df['shape_score'].dropna()
    logger.info(f"Shape score statistics ({method}):")
    logger.info(f"  Mean: {scores.mean():.4f}")
    logger.info(f"  Median: {scores.median():.4f}")
    logger.info(f"  Min: {scores.min():.4f}")
    logger.info(f"  Max: {scores.max():.4f}")


@app.command("batch")
def batch_command(
    cluster_file: Path = typer.Option(
        ..., "-c", "--clusters",
        help="Input cluster file"
    ),
    tss_file: Path = typer.Option(
        ..., "-t", "--tss",
        help="Input TSS table"
    ),
    output_prefix: str = typer.Option(
        ..., "-o", "--output-prefix",
        help="Output file prefix"
    ),
    method: str = typer.Option(
        "PSS", "-m", "--method",
        help="Shape calculation method: PSS or SI"
    ),
    processes: int = typer.Option(
        1, "-p", "--processes",
        help="Number of processes"
    ),
):
    """
    Calculate shape scores for all samples in a TSS table.

    Creates one output file per sample: {prefix}.{sample}.shape.tsv

    Example:
        tsspy shapeCluster batch -c clusters.tsv -t tss.tsv -o shape -m PSS
    """
    logger.info(f"Loading cluster file: {cluster_file}")
    cluster_df = pd.read_csv(cluster_file, sep='\t')

    logger.info(f"Loading TSS file: {tss_file}")
    tss_df = pd.read_csv(tss_file, sep='\t')

    sample_cols = [c for c in tss_df.columns if c not in ['chr', 'pos', 'strand']]
    logger.info(f"Samples: {sample_cols}")

    for sample in sample_cols:
        logger.info(f"Processing sample: {sample}")

        result_df = calculate_shape_scores(
            cluster_df, tss_df, sample, method=method, n_processes=processes
        )

        output_file = f"{output_prefix}.{sample}.shape.tsv"
        result_df.to_csv(output_file, sep='\t', index=False)
        logger.info(f"Saved: {output_file}")


@app.command("classify")
def classify_command(
    shape_file: Path = typer.Option(
        ..., "-i", "--input",
        help="Input shape score file"
    ),
    output_file: Path = typer.Option(
        ..., "-o", "--output",
        help="Output classified file"
    ),
    method: str = typer.Option(
        "PSS", "-m", "--method",
        help="Method used for shape scores (PSS or SI)"
    ),
    threshold: float = typer.Option(
        None, "-t", "--threshold",
        help="Custom threshold for sharp/broad classification"
    ),
):
    """
    Classify promoters as sharp or broad based on shape scores.

    Default thresholds:
        - PSS: < 5 = sharp, >= 5 = broad
        - SI: > 1.5 = sharp, <= 1.5 = broad

    Example:
        tsspy shapeCluster classify -i shape.tsv -o classified.tsv -m PSS
    """
    logger.info(f"Loading shape file: {shape_file}")
    df = pd.read_csv(shape_file, sep='\t')

    if 'shape_score' not in df.columns:
        raise typer.BadParameter("Input file must have 'shape_score' column")

    # Set default thresholds
    if threshold is None:
        if method.upper() == 'PSS':
            threshold = 5.0
        else:  # SI
            threshold = 1.5

    # Classify promoters
    if method.upper() == 'PSS':
        # Lower PSS = sharper
        df['promoter_class'] = df['shape_score'].apply(
            lambda x: 'sharp' if x < threshold else 'broad'
        )
    else:  # SI
        # Higher SI = sharper
        df['promoter_class'] = df['shape_score'].apply(
            lambda x: 'sharp' if x > threshold else 'broad'
        )

    # Save results
    df.to_csv(output_file, sep='\t', index=False)
    logger.info(f"Saved classified results to {output_file}")

    # Report statistics
    class_counts = df['promoter_class'].value_counts()
    logger.info(f"Classification results (threshold={threshold}):")
    for cls, count in class_counts.items():
        pct = count / len(df) * 100
        logger.info(f"  {cls}: {count} ({pct:.1f}%)")


if __name__ == '__main__':
    app()
