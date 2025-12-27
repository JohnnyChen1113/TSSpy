#!/usr/bin/env python3
"""
Consensus Cluster - Aggregate tag clusters across samples into consensus clusters
"""

import typer
import pandas as pd
import numpy as np
from typing import Optional, List, Dict
from pathlib import Path
import logging
from collections import defaultdict

logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

__version__ = "0.1.0"

app = typer.Typer(help=f"Consensus clustering across samples (v{__version__})")


def load_tag_clusters(cluster_files: List[str], sample_names: List[str]) -> Dict[str, pd.DataFrame]:
    """
    Load tag cluster files for multiple samples.

    Args:
        cluster_files: List of paths to cluster TSV files
        sample_names: Corresponding sample names

    Returns:
        Dictionary mapping sample names to cluster DataFrames
    """
    clusters = {}
    for f, name in zip(cluster_files, sample_names):
        df = pd.read_csv(f, sep='\t')
        # Ensure required columns exist
        required_cols = ['cluster', 'chr', 'start', 'end', 'strand', 'dominant_tss']
        missing = [c for c in required_cols if c not in df.columns]
        if missing:
            logger.warning(f"Sample {name}: missing columns {missing}")
            continue
        clusters[name] = df
        logger.info(f"Loaded {len(df)} clusters for sample {name}")
    return clusters


def find_overlapping_clusters(clusters_dict: Dict[str, pd.DataFrame],
                               dominant_distance: int = 50) -> pd.DataFrame:
    """
    Find consensus clusters by grouping clusters with nearby dominant TSS.

    Clusters from different samples are considered to belong to the same
    consensus cluster if their dominant TSS positions are within the
    specified distance.

    Args:
        clusters_dict: Dictionary of sample -> cluster DataFrame
        dominant_distance: Maximum distance between dominant TSS positions

    Returns:
        DataFrame with consensus cluster information
    """
    # Collect all clusters with their sample information
    all_clusters = []
    for sample, df in clusters_dict.items():
        df_copy = df.copy()
        df_copy['sample'] = sample
        all_clusters.append(df_copy)

    if not all_clusters:
        return pd.DataFrame()

    combined = pd.concat(all_clusters, ignore_index=True)

    # Sort by chromosome, strand, and dominant TSS position
    combined = combined.sort_values(['chr', 'strand', 'dominant_tss']).reset_index(drop=True)

    # Assign consensus cluster IDs
    consensus_id = 0
    consensus_ids = []
    prev_chr = None
    prev_strand = None
    prev_dom_tss = None

    for idx, row in combined.iterrows():
        chr_name = row['chr']
        strand = row['strand']
        dom_tss = row['dominant_tss']

        # Check if this cluster belongs to a new consensus group
        if (prev_chr != chr_name or
            prev_strand != strand or
            prev_dom_tss is None or
            abs(dom_tss - prev_dom_tss) > dominant_distance):
            consensus_id += 1

        consensus_ids.append(consensus_id)
        prev_chr = chr_name
        prev_strand = strand
        prev_dom_tss = dom_tss

    combined['consensus_cluster'] = consensus_ids

    return combined


def aggregate_consensus_clusters(combined: pd.DataFrame,
                                  sample_names: List[str]) -> pd.DataFrame:
    """
    Aggregate consensus cluster information across samples.

    For each consensus cluster, compute:
    - Combined boundaries (union of all sample boundaries)
    - Dominant TSS (most common or highest signal)
    - Signal statistics per sample

    Args:
        combined: DataFrame with consensus_cluster column
        sample_names: List of sample names

    Returns:
        Aggregated consensus cluster DataFrame
    """
    results = []

    for cc_id, group in combined.groupby('consensus_cluster'):
        # Basic info (from first row, they should be similar)
        chr_name = group['chr'].iloc[0]
        strand = group['strand'].iloc[0]

        # Combined boundaries
        start = group['start'].min()
        end = group['end'].max()

        # Dominant TSS: use the position with highest total tags
        if 'tags' in group.columns:
            dom_idx = group['tags'].idxmax()
            dominant_tss = group.loc[dom_idx, 'dominant_tss']
        else:
            # Use the median dominant TSS position
            dominant_tss = int(group['dominant_tss'].median())

        # Collect per-sample statistics
        row = {
            'consensus_cluster': cc_id,
            'chr': chr_name,
            'start': start,
            'end': end,
            'strand': strand,
            'dominant_tss': dominant_tss,
            'width': end - start + 1,
            'n_samples': group['sample'].nunique(),
        }

        # Add sample-specific info
        for sample in sample_names:
            sample_data = group[group['sample'] == sample]
            if len(sample_data) > 0:
                row[f'{sample}.tags'] = sample_data['tags'].sum() if 'tags' in sample_data.columns else 0
                row[f'{sample}.dominant_tss'] = sample_data['dominant_tss'].iloc[0]
            else:
                row[f'{sample}.tags'] = 0
                row[f'{sample}.dominant_tss'] = np.nan

        results.append(row)

    return pd.DataFrame(results)


def compute_consensus_from_tss_table(tss_df: pd.DataFrame,
                                      clustered_dfs: Dict[str, pd.DataFrame],
                                      dominant_distance: int = 50) -> pd.DataFrame:
    """
    Alternative approach: create consensus clusters directly from TSS table
    by finding regions where multiple samples have clusters.

    Args:
        tss_df: Raw TSS table
        clustered_dfs: Per-sample clustered data
        dominant_distance: Distance threshold for merging

    Returns:
        Consensus cluster DataFrame with per-sample tags
    """
    # Find consensus clusters
    combined = find_overlapping_clusters(clustered_dfs, dominant_distance)

    if combined.empty:
        return pd.DataFrame()

    sample_names = list(clustered_dfs.keys())
    consensus_df = aggregate_consensus_clusters(combined, sample_names)

    return consensus_df


@app.command("cluster")
def consensus_cluster_command(
    input_files: List[Path] = typer.Option(
        ..., "-i", "--input",
        help="Input tag cluster files (one per sample)"
    ),
    sample_names: str = typer.Option(
        ..., "-n", "--sample-names",
        help="Sample names (space-separated)"
    ),
    output_file: Path = typer.Option(
        ..., "-o", "--output",
        help="Output consensus cluster file"
    ),
    distance: int = typer.Option(
        50, "-d", "--distance",
        help="Maximum distance between dominant TSS positions to merge clusters"
    ),
):
    """
    Create consensus clusters from per-sample tag clusters.

    Example:
        tsspy consensusCluster cluster \\
            -i control.clusters.tsv -i treat.clusters.tsv \\
            -n "control treat" -o consensus.tsv -d 50
    """
    names_list = sample_names.split()
    files_list = [str(f) for f in input_files]

    if len(names_list) != len(files_list):
        raise typer.BadParameter(
            f"Number of sample names ({len(names_list)}) must match number of input files ({len(files_list)})"
        )

    logger.info(f"Loading {len(files_list)} cluster files")

    # Load cluster files
    clusters_dict = load_tag_clusters(files_list, names_list)

    if not clusters_dict:
        logger.error("No valid cluster files loaded")
        raise typer.Exit(1)

    # Find consensus clusters
    logger.info(f"Finding consensus clusters (distance threshold: {distance} bp)")
    combined = find_overlapping_clusters(clusters_dict, distance)

    if combined.empty:
        logger.error("No clusters found")
        raise typer.Exit(1)

    # Aggregate
    consensus_df = aggregate_consensus_clusters(combined, names_list)

    # Save
    consensus_df.to_csv(output_file, sep='\t', index=False)
    logger.info(f"Saved {len(consensus_df)} consensus clusters to {output_file}")


@app.command("from-tss")
def consensus_from_tss_command(
    tss_input: Path = typer.Option(
        ..., "-t", "--tss-input",
        help="Input TSS table (merged/normalized)"
    ),
    output_file: Path = typer.Option(
        ..., "-o", "--output",
        help="Output consensus cluster file"
    ),
    peak_distance: int = typer.Option(
        100, "--peak-distance",
        help="Minimum distance between peaks for clustering"
    ),
    extension_distance: int = typer.Option(
        30, "--extension-distance",
        help="Cluster boundary extension distance"
    ),
    local_threshold: float = typer.Option(
        0.02, "--local-threshold",
        help="Local filtering threshold (fraction of peak TPM)"
    ),
    cluster_threshold: float = typer.Option(
        1.0, "--cluster-threshold",
        help="Minimum cluster TPM threshold"
    ),
    consensus_distance: int = typer.Option(
        50, "--consensus-distance",
        help="Distance threshold for consensus clustering"
    ),
):
    """
    Cluster each sample and then create consensus clusters.

    This performs per-sample clustering followed by consensus aggregation.

    Example:
        tsspy consensusCluster from-tss -t normalized.tsv -o consensus.tsv
    """
    from TSSpy.clustering import cluster_by_peak

    logger.info(f"Reading TSS table: {tss_input}")
    df = pd.read_csv(tss_input, sep='\t')

    sample_cols = [c for c in df.columns if c not in ['chr', 'pos', 'strand']]
    logger.info(f"Samples: {sample_cols}")

    # Cluster each sample
    all_clusters = {}

    for sample in sample_cols:
        logger.info(f"Clustering sample: {sample}")

        # Calculate TPM for this sample
        total = df[sample].sum()
        if total == 0:
            logger.warning(f"Sample {sample} has zero total counts, skipping")
            continue

        sample_df = df[['chr', 'pos', 'strand', sample]].copy()
        sample_df['TPM'] = sample_df[sample] / total * 1e6
        sample_df = sample_df[sample_df['TPM'] > 0]

        # Cluster by chromosome and strand
        sample_clusters = []
        for (chr_name, strand), group in sample_df.groupby(['chr', 'strand']):
            group = group.sort_values('pos').reset_index(drop=True)
            clusters = cluster_by_peak(
                group, peak_distance, local_threshold,
                extension_distance, sample, cluster_threshold
            )
            if clusters is not None and not clusters.empty:
                clusters['chr'] = chr_name
                clusters['strand'] = strand
                sample_clusters.append(clusters)

        if sample_clusters:
            all_clusters[sample] = pd.concat(sample_clusters, ignore_index=True)
            logger.info(f"Sample {sample}: {len(all_clusters[sample])} clusters")

    if not all_clusters:
        logger.error("No clusters found for any sample")
        raise typer.Exit(1)

    # Create consensus clusters
    logger.info(f"Creating consensus clusters (distance: {consensus_distance} bp)")
    combined = find_overlapping_clusters(all_clusters, consensus_distance)
    consensus_df = aggregate_consensus_clusters(combined, sample_cols)

    # Save
    consensus_df.to_csv(output_file, sep='\t', index=False)
    logger.info(f"Saved {len(consensus_df)} consensus clusters to {output_file}")


if __name__ == '__main__':
    app()
