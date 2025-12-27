#!/usr/bin/env python3
"""
Gene Assign - Assign TSS clusters to genes based on genomic annotations

Supports GTF and GFF annotation formats.
"""

import typer
import pandas as pd
import numpy as np
from typing import Optional, List, Dict, Tuple
from pathlib import Path
import logging
import re

logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

__version__ = "0.1.0"

app = typer.Typer(help=f"Assign TSS clusters to genes (v{__version__})")


def parse_gtf_attributes(attr_string: str) -> Dict[str, str]:
    """Parse GTF attribute string into dictionary."""
    attrs = {}
    # GTF format: key "value"; key "value";
    pattern = r'(\w+)\s+"([^"]+)"'
    for match in re.finditer(pattern, attr_string):
        attrs[match.group(1)] = match.group(2)
    return attrs


def parse_gff_attributes(attr_string: str) -> Dict[str, str]:
    """Parse GFF3 attribute string into dictionary."""
    attrs = {}
    # GFF3 format: key=value;key=value
    for item in attr_string.split(';'):
        item = item.strip()
        if '=' in item:
            key, value = item.split('=', 1)
            # Handle URL encoding
            value = value.replace('%2C', ',').replace('%3B', ';')
            attrs[key] = value
    return attrs


def load_annotation(annotation_file: str, annotation_type: str = 'genes') -> pd.DataFrame:
    """
    Load gene annotation from GTF or GFF file.

    Args:
        annotation_file: Path to GTF/GFF file
        annotation_type: 'genes' for CDS-based, 'transcript' for transcript-based

    Returns:
        DataFrame with gene annotation (chr, start, end, strand, gene_id, gene_name)
    """
    logger.info(f"Loading annotation from: {annotation_file}")

    # Detect file format
    is_gtf = annotation_file.lower().endswith('.gtf')

    records = []

    with open(annotation_file, 'r') as f:
        for line in f:
            if line.startswith('#'):
                continue

            fields = line.strip().split('\t')
            if len(fields) < 9:
                continue

            chr_name = fields[0]
            feature_type = fields[2]
            start = int(fields[3])  # 1-based
            end = int(fields[4])
            strand = fields[6]
            attrs_str = fields[8]

            # Filter by feature type
            if annotation_type == 'genes':
                # Use CDS or gene features
                if feature_type not in ['CDS', 'gene', 'mRNA']:
                    continue
            else:  # transcript
                if feature_type not in ['transcript', 'mRNA', 'gene']:
                    continue

            # Parse attributes
            if is_gtf:
                attrs = parse_gtf_attributes(attrs_str)
                gene_id = attrs.get('gene_id', attrs.get('gene_name', ''))
                gene_name = attrs.get('gene_name', gene_id)
            else:
                attrs = parse_gff_attributes(attrs_str)
                gene_id = attrs.get('ID', attrs.get('Name', attrs.get('gene', '')))
                gene_name = attrs.get('Name', attrs.get('gene', gene_id))
                # Handle Parent attribute for mRNA/CDS features
                if 'Parent' in attrs and not gene_id:
                    gene_id = attrs['Parent']

            if gene_id:
                records.append({
                    'chr': chr_name,
                    'start': start,
                    'end': end,
                    'strand': strand,
                    'gene_id': gene_id,
                    'gene_name': gene_name,
                    'feature_type': feature_type
                })

    df = pd.DataFrame(records)

    # Aggregate by gene (get min start and max end for each gene)
    if len(df) > 0:
        gene_df = df.groupby(['chr', 'strand', 'gene_id', 'gene_name']).agg({
            'start': 'min',
            'end': 'max'
        }).reset_index()

        # For plus strand genes, the 5' end is the start
        # For minus strand genes, the 5' end is the end
        gene_df['five_prime'] = np.where(
            gene_df['strand'] == '+',
            gene_df['start'],
            gene_df['end']
        )

        logger.info(f"Loaded {len(gene_df)} genes")
        return gene_df

    return pd.DataFrame()


def assign_cluster_to_gene(cluster: pd.Series,
                            gene_df: pd.DataFrame,
                            upstream: int = 1000,
                            downstream: int = 0,
                            upstream_overlap: int = 500) -> Tuple[str, bool]:
    """
    Assign a single cluster to its downstream gene.

    Assignment rules (from TSSr):
    1. The dominant TSS must be upstream of the gene's 5' end (within upstream distance)
    2. If the TSS overlaps with an upstream gene's coding region, it must be within
       upstream_overlap distance of that gene's 3' end

    Args:
        cluster: Series with cluster information (must have chr, strand, dominant_tss)
        gene_df: Gene annotation DataFrame
        upstream: Maximum distance upstream of gene 5' end
        downstream: Maximum distance downstream of gene 5' end (into gene body)
        upstream_overlap: Max overlap with upstream gene's coding region

    Returns:
        Tuple of (gene_id, is_in_coding_region)
    """
    chr_name = cluster['chr']
    strand = cluster['strand']
    dom_tss = cluster['dominant_tss']

    # Filter genes on same chromosome and strand
    genes = gene_df[(gene_df['chr'] == chr_name) & (gene_df['strand'] == strand)]

    if genes.empty:
        return '', False

    best_gene = ''
    in_coding = False
    best_distance = float('inf')

    for _, gene in genes.iterrows():
        gene_id = gene['gene_id']
        gene_start = gene['start']
        gene_end = gene['end']
        five_prime = gene['five_prime']

        if strand == '+':
            # Plus strand: TSS should be upstream (smaller position) of gene start
            distance = five_prime - dom_tss

            # Check if TSS is in valid range
            if -downstream <= distance <= upstream:
                # Check if TSS overlaps with coding region
                is_inside = gene_start <= dom_tss <= gene_end

                if is_inside:
                    # TSS is inside the gene - this might be internal TSS
                    in_coding = True

                if distance >= 0 and distance < best_distance:
                    best_gene = gene_id
                    best_distance = distance
                elif distance < 0 and abs(distance) <= downstream and not best_gene:
                    # TSS is slightly downstream (inside gene)
                    best_gene = gene_id
                    in_coding = True

        else:
            # Minus strand: TSS should be upstream (larger position) of gene end
            distance = dom_tss - five_prime

            if -downstream <= distance <= upstream:
                is_inside = gene_start <= dom_tss <= gene_end

                if is_inside:
                    in_coding = True

                if distance >= 0 and distance < best_distance:
                    best_gene = gene_id
                    best_distance = distance
                elif distance < 0 and abs(distance) <= downstream and not best_gene:
                    best_gene = gene_id
                    in_coding = True

    return best_gene, in_coding


def assign_clusters_to_genes(cluster_df: pd.DataFrame,
                              gene_df: pd.DataFrame,
                              upstream: int = 1000,
                              downstream: int = 0,
                              upstream_overlap: int = 500,
                              filter_threshold: float = 0.0) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """
    Assign all clusters to genes.

    Args:
        cluster_df: Cluster DataFrame
        gene_df: Gene annotation DataFrame
        upstream: Max upstream distance
        downstream: Max downstream distance
        upstream_overlap: Max overlap with upstream gene
        filter_threshold: Minimum cluster TPM to include

    Returns:
        Tuple of (assigned_clusters, unassigned_clusters)
    """
    # Filter clusters if threshold specified
    if filter_threshold > 0 and 'tags_TPM' in cluster_df.columns:
        cluster_df = cluster_df[cluster_df['tags_TPM'] >= filter_threshold].copy()
    elif filter_threshold > 0 and 'tags' in cluster_df.columns:
        cluster_df = cluster_df[cluster_df['tags'] >= filter_threshold].copy()

    assignments = []
    for idx, cluster in cluster_df.iterrows():
        gene_id, in_coding = assign_cluster_to_gene(
            cluster, gene_df, upstream, downstream, upstream_overlap
        )
        assignments.append({
            'gene': gene_id,
            'in_coding': in_coding
        })

    assign_df = pd.DataFrame(assignments)
    result = pd.concat([cluster_df.reset_index(drop=True), assign_df], axis=1)

    # Split into assigned and unassigned
    assigned = result[result['gene'] != ''].copy()
    unassigned = result[result['gene'] == ''].copy()

    logger.info(f"Assigned: {len(assigned)} clusters, Unassigned: {len(unassigned)} clusters")

    return assigned, unassigned


@app.command("assign")
def assign_command(
    cluster_file: Path = typer.Option(
        ..., "-c", "--clusters",
        help="Input cluster file"
    ),
    annotation_file: Path = typer.Option(
        ..., "-a", "--annotation",
        help="Gene annotation file (GTF or GFF format)"
    ),
    output_assigned: Path = typer.Option(
        ..., "-o", "--output",
        help="Output file for assigned clusters"
    ),
    output_unassigned: Optional[Path] = typer.Option(
        None, "-u", "--unassigned",
        help="Output file for unassigned clusters (optional)"
    ),
    annotation_type: str = typer.Option(
        "genes", "--annotation-type",
        help="Annotation type: 'genes' (CDS-based) or 'transcript'"
    ),
    upstream: int = typer.Option(
        1000, "--upstream",
        help="Maximum distance upstream of gene 5' end"
    ),
    downstream: int = typer.Option(
        0, "--downstream",
        help="Maximum distance downstream (into gene body)"
    ),
    filter_threshold: float = typer.Option(
        0.02, "--filter-threshold",
        help="Minimum cluster TPM threshold"
    ),
):
    """
    Assign TSS clusters to downstream genes.

    A cluster is assigned to a gene if its dominant TSS is within the
    specified upstream distance of the gene's transcription start site.

    Example:
        tsspy geneAssign assign -c clusters.tsv -a annotation.gtf -o assigned.tsv
    """
    logger.info(f"Loading clusters: {cluster_file}")
    cluster_df = pd.read_csv(cluster_file, sep='\t')

    # Load annotation
    gene_df = load_annotation(str(annotation_file), annotation_type)

    if gene_df.empty:
        logger.error("No genes found in annotation file")
        raise typer.Exit(1)

    logger.info(f"Assigning clusters (upstream={upstream}, downstream={downstream})")

    assigned, unassigned = assign_clusters_to_genes(
        cluster_df, gene_df,
        upstream=upstream,
        downstream=downstream,
        filter_threshold=filter_threshold
    )

    # Save assigned clusters
    assigned.to_csv(output_assigned, sep='\t', index=False)
    logger.info(f"Saved assigned clusters to {output_assigned}")

    # Save unassigned clusters if requested
    if output_unassigned:
        unassigned.to_csv(output_unassigned, sep='\t', index=False)
        logger.info(f"Saved unassigned clusters to {output_unassigned}")

    # Report statistics
    if 'gene' in assigned.columns:
        unique_genes = assigned['gene'].nunique()
        logger.info(f"Total assigned clusters: {len(assigned)}")
        logger.info(f"Unique genes with assigned clusters: {unique_genes}")

        if 'in_coding' in assigned.columns:
            in_coding_count = assigned['in_coding'].sum()
            logger.info(f"Clusters inside coding regions: {in_coding_count}")


@app.command("summary")
def summary_command(
    assigned_file: Path = typer.Option(
        ..., "-i", "--input",
        help="Input assigned clusters file"
    ),
    output_file: Path = typer.Option(
        ..., "-o", "--output",
        help="Output gene-level summary file"
    ),
):
    """
    Create a gene-level summary from assigned clusters.

    Aggregates cluster information by gene, reporting:
    - Number of clusters per gene
    - Total and max tags
    - Primary cluster (highest tags)

    Example:
        tsspy geneAssign summary -i assigned.tsv -o gene_summary.tsv
    """
    logger.info(f"Loading assigned clusters: {assigned_file}")
    df = pd.read_csv(assigned_file, sep='\t')

    if 'gene' not in df.columns:
        raise typer.BadParameter("Input file must have 'gene' column")

    # Determine tag column
    tag_col = 'tags_TPM' if 'tags_TPM' in df.columns else 'tags'

    # Group by gene
    summary = df.groupby('gene').agg({
        'chr': 'first',
        'strand': 'first',
        'dominant_tss': lambda x: list(x),
        tag_col: ['sum', 'max', 'count']
    }).reset_index()

    # Flatten column names
    summary.columns = ['gene', 'chr', 'strand', 'dominant_tss_positions',
                       'total_tags', 'max_cluster_tags', 'n_clusters']

    # Convert list to string for output
    summary['dominant_tss_positions'] = summary['dominant_tss_positions'].apply(
        lambda x: ','.join(map(str, x))
    )

    # Sort by total tags
    summary = summary.sort_values('total_tags', ascending=False)

    # Save
    summary.to_csv(output_file, sep='\t', index=False)
    logger.info(f"Saved gene summary to {output_file}")
    logger.info(f"Total genes: {len(summary)}")


@app.command("filter")
def filter_command(
    assigned_file: Path = typer.Option(
        ..., "-i", "--input",
        help="Input assigned clusters file"
    ),
    output_file: Path = typer.Option(
        ..., "-o", "--output",
        help="Output filtered file"
    ),
    keep_primary: bool = typer.Option(
        False, "--keep-primary",
        help="Keep only the primary (highest tags) cluster per gene"
    ),
    exclude_internal: bool = typer.Option(
        False, "--exclude-internal",
        help="Exclude clusters inside coding regions"
    ),
    min_tags: float = typer.Option(
        0.0, "--min-tags",
        help="Minimum tags threshold"
    ),
):
    """
    Filter assigned clusters based on various criteria.

    Example:
        tsspy geneAssign filter -i assigned.tsv -o filtered.tsv --keep-primary
    """
    logger.info(f"Loading assigned clusters: {assigned_file}")
    df = pd.read_csv(assigned_file, sep='\t')

    original_count = len(df)

    # Filter by internal status
    if exclude_internal and 'in_coding' in df.columns:
        df = df[~df['in_coding']]
        logger.info(f"Excluded internal clusters: {original_count} -> {len(df)}")

    # Filter by min tags
    tag_col = 'tags_TPM' if 'tags_TPM' in df.columns else 'tags'
    if min_tags > 0:
        df = df[df[tag_col] >= min_tags]
        logger.info(f"Filtered by min tags ({min_tags}): {len(df)} remaining")

    # Keep only primary cluster per gene
    if keep_primary and 'gene' in df.columns:
        df = df.loc[df.groupby('gene')[tag_col].idxmax()]
        logger.info(f"Kept primary clusters only: {len(df)} remaining")

    # Save
    df.to_csv(output_file, sep='\t', index=False)
    logger.info(f"Saved filtered clusters to {output_file}")


if __name__ == '__main__':
    app()
