#!/usr/bin/env python3
"""
Assign TSS clusters to genes — strict parity with TSSr 0.99.6 annotateCluster().

Mirrors:
  TSSr::annotateCluster()       (AnnotationMethods.R)
  TSSr::.assign2gene()          (AnnotationFunctions.R)

Algorithm:
  1. Build promoter regions per gene with neighbor-aware up/down distances.
  2. Overlap each cluster's dominant_tss with promoter regions → assign `gene`
     (deduplicated by first hit per cluster).
  3. Overlap each cluster's dominant_tss with raw gene bodies → assign `inCoding`.
  4. If filterCluster: within each (gene-or-inCoding) group, drop clusters that
     are downstream of the dominant cluster AND below `filterClusterThreshold *
     max(tags)`.

Outputs:
  - assigned    : clusters with non-NA gene
  - unassigned  : clusters with NA gene
  - filtered    : assigned+inCoding rows after the filter, plus clusters with
                  no gene+no inCoding (unfiltered)
"""

from __future__ import annotations
import typer
import pandas as pd
import numpy as np
from pathlib import Path
from typing import List, Optional, Dict, Tuple
from urllib.parse import unquote
import logging

logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

__version__ = "0.18.0"
app = typer.Typer(help=f"Assign TSS clusters to genes (v{__version__})")


# -----------------------------------------------------------------------------
# GFF loading (mirrors GenomicFeatures::genes() output)
# -----------------------------------------------------------------------------

# Gene-like GFF feature types accepted as "gene" (matches what TSSr's
# makeTxDbFromGFF + genes(txdb) treats as gene features in SGD-style GFFs).
GENE_FEATURE_TYPES = {
    'gene', 'tRNA_gene', 'transposable_element_gene', 'snoRNA_gene',
    'rRNA_gene', 'ncRNA_gene', 'pseudogene', 'snRNA_gene',
    'telomerase_RNA_gene',
}


def load_genes_from_gff(gff_path: str) -> pd.DataFrame:
    """
    Read a GFF3 file and return one row per gene-like feature with columns:
    seqnames, start, end, strand, width, gene_id.
    Matches the structure of TSSr's `as.data.frame(genes(txdb))`.

    Accepts the full set of SGD-style gene-like feature types
    (gene, tRNA_gene, snoRNA_gene, etc.) and URL-decodes the ID attribute
    (e.g. `tP%28UGG%29A` -> `tP(UGG)A`).
    """
    rows = []
    with open(gff_path) as fh:
        for line in fh:
            if not line or line.startswith('#'):
                continue
            parts = line.rstrip('\n').split('\t')
            if len(parts) < 9:
                continue
            if parts[2] not in GENE_FEATURE_TYPES:
                continue
            chrom, start, end, strand, attrs = parts[0], int(parts[3]), int(parts[4]), parts[6], parts[8]
            gene_id = ''
            for kv in attrs.split(';'):
                kv = kv.strip()
                if kv.startswith('ID='):
                    gene_id = unquote(kv[3:])
                    break
            if not gene_id:
                continue
            rows.append((chrom, start, end, strand, end - start + 1, gene_id))
    df = pd.DataFrame(rows, columns=['seqnames', 'start', 'end', 'strand', 'width', 'gene_id'])
    return df


# -----------------------------------------------------------------------------
# Promoter region construction (TSSr .assign2gene up/down logic)
# -----------------------------------------------------------------------------

def _promoter_regions(ref_sub: pd.DataFrame, strand: str,
                      upstream: int, upstream_overlap: int, downstream: int) -> pd.DataFrame:
    """
    Per (chr, strand): construct promoter regions with neighbor-aware up/down
    distances. Input must be a single-(chr,strand) slice with columns
    seqnames, start, end, strand, width, gene_id.
    Returns a copy with start/end overwritten to be the promoter range.
    """
    df = ref_sub.copy()
    if strand == '+':
        # R's data.table::setorder is stable; pandas default quicksort isn't.
        # Stable sort is required to match TSSr when multiple genes share a start.
        df = df.sort_values('start', kind='stable').reset_index(drop=True)
        end_b = df['end'].shift(1, fill_value=0).to_numpy()       # prev gene's end
        width = df['width'].shift(1, fill_value=1000).to_numpy()  # prev gene's width
        dis = df['start'].to_numpy() - end_b
        up = np.where(dis > upstream, upstream,
              np.where(dis + width <= upstream_overlap, dis + width - 1,
               np.where(dis < upstream_overlap, upstream_overlap, dis)))
        start_a = df['start'].shift(-1, fill_value=1000).to_numpy()  # next gene's start
        dis_start = start_a - df['start'].to_numpy()
        down_up = np.roll(up, -1)  # next gene's up
        down_up[-1] = 1000
        down = np.where(dis_start > downstream + down_up, downstream,
                np.where(dis_start < down_up, 0,
                         dis_start - down_up))
        new_end = df['start'].to_numpy() + down
        new_start = new_end - down - up + 1
    else:  # minus strand
        df = df.sort_values('end', kind='stable').reset_index(drop=True)
        # end.b = lead(start, 1, fill=0); last row overrides to end + 1000
        end_b = df['start'].shift(-1, fill_value=0).to_numpy()
        end_b[-1] = df['end'].iloc[-1] + 1000
        width = df['width'].shift(-1, fill_value=1000).to_numpy()
        dis = end_b - df['end'].to_numpy()
        up = np.where(dis > upstream, upstream,
              np.where(dis + width <= upstream_overlap, dis + width - 1,
               np.where(dis < upstream_overlap, upstream_overlap, dis)))
        start_a = df['end'].shift(1, fill_value=1000).to_numpy()  # prev gene's end (in this sort order)
        dis_start = df['end'].to_numpy() - start_a
        down_up = np.roll(up, 1)  # prev gene's up
        down_up[0] = 1000
        down = np.where(dis_start >= downstream + down_up, downstream,
                np.where(dis_start < down_up, 0,
                         dis_start - down_up))
        new_start = df['end'].to_numpy() - down
        new_end = new_start + down + up - 1
    df['start'] = new_start.astype(np.int64)
    df['end'] = new_end.astype(np.int64)
    return df


# -----------------------------------------------------------------------------
# Overlap helpers (single-point queries against ranges)
# -----------------------------------------------------------------------------

def _first_overlapping_gene(positions: np.ndarray, ranges: pd.DataFrame) -> List[Optional[str]]:
    """
    For each position (single-base), find the first range (by row order) that
    contains it inclusively. Returns the gene_id of that range, or None.
    Mirrors TSSr's `findOverlaps` then dedup-by-queryHits behaviour.
    """
    starts = ranges['start'].to_numpy()
    ends = ranges['end'].to_numpy()
    gene_ids = ranges['gene_id'].to_numpy()
    out: List[Optional[str]] = []
    for p in positions:
        # Naive O(N) for clarity; cluster counts are in thousands so fine
        hits = np.where((starts <= p) & (ends >= p))[0]
        if len(hits) == 0:
            out.append(None)
        else:
            out.append(str(gene_ids[hits[0]]))
    return out


# -----------------------------------------------------------------------------
# Driver: assign one sample's clusters to genes
# -----------------------------------------------------------------------------

def assign_one_sample(cluster_df: pd.DataFrame,
                      ref_genes: pd.DataFrame,
                      upstream: int = 1000,
                      upstream_overlap: int = 500,
                      downstream: int = 0,
                      filter_cluster: bool = True) -> pd.DataFrame:
    """
    For one sample's clusters, assign gene + (optional) inCoding, sort by cluster.
    Output mirrors TSSr's per-sample asn DataFrame (before split into assigned/
    unassigned/filtered).
    """
    # Stratify both inputs by (chr, strand) — only matching strata interact.
    cs_chr_strand = list(cluster_df.groupby(['chr', 'strand'], sort=False).groups.keys())
    ref_chr_strand = set(map(tuple, ref_genes[['seqnames', 'strand']].drop_duplicates().to_numpy()))

    rows = []
    for (chrom, strand), idx in cluster_df.groupby(['chr', 'strand'], sort=False).groups.items():
        cs_grp = cluster_df.loc[idx].copy()
        if (chrom, strand) not in ref_chr_strand:
            cs_grp['gene'] = None
            if filter_cluster:
                cs_grp['inCoding'] = None
            rows.append(cs_grp)
            continue

        ref_sub = ref_genes[(ref_genes['seqnames'] == chrom) & (ref_genes['strand'] == strand)].copy()

        # Promoter overlap
        prom = _promoter_regions(ref_sub, strand, upstream, upstream_overlap, downstream)
        cs_grp['gene'] = _first_overlapping_gene(cs_grp['dominant_tss'].to_numpy(),
                                                  prom[['start', 'end', 'gene_id']])
        # Coding-body overlap (raw gene start/end, no promoter expansion)
        if filter_cluster:
            ref_coding = ref_sub.sort_values('start', kind='stable').reset_index(drop=True)
            cs_grp['inCoding'] = _first_overlapping_gene(cs_grp['dominant_tss'].to_numpy(),
                                                          ref_coding[['start', 'end', 'gene_id']])
        rows.append(cs_grp)

    if not rows:
        return cluster_df.copy()
    out = pd.concat(rows, ignore_index=False)
    out = out.sort_values('cluster').reset_index(drop=True)
    return out


def filter_assigned(asn: pd.DataFrame, filter_cluster_threshold: float) -> pd.DataFrame:
    """
    TSSr's filter step: within each (gene-or-inCoding) group, drop clusters
    that are downstream of the dominant cluster AND below threshold.
    Clusters with both gene and inCoding NA pass through unfiltered.
    """
    m = asn[asn['gene'].isna() & asn['inCoding'].isna()].copy()
    n = asn[asn['gene'].notna() | asn['inCoding'].notna()].copy()
    if n.empty:
        return m

    # If both gene and inCoding set, drop inCoding (gene assignment wins)
    n['inCoding'] = np.where(n['gene'].notna() & n['inCoding'].notna(), None, n['inCoding'])
    n['r'] = n['gene'].where(n['gene'].notna(), n['inCoding'])

    kept = []
    for r_val, grp in n.groupby('r', sort=False):
        max_idx = grp['tags'].idxmax()
        max_tags = float(grp.at[max_idx, 'tags'])
        dom_of_max = int(grp.at[max_idx, 'dominant_tss'])
        thr = max_tags * filter_cluster_threshold
        if grp['strand'].iloc[0] == '+':
            f = ~((grp['dominant_tss'] > dom_of_max) & (grp['tags'] < thr))
        else:
            f = ~((grp['dominant_tss'] < dom_of_max) & (grp['tags'] < thr))
        kept.append(grp[f])
    n_filtered = pd.concat(kept, ignore_index=True) if kept else pd.DataFrame(columns=n.columns)
    n_filtered = n_filtered.drop(columns=['r'], errors='ignore')

    # TSSr's `rbind(m[,seq(12)], new[,seq(12)])` takes the first 12 columns of each.
    keep_cols = [c for c in m.columns if c != 'r']
    out = pd.concat([m[keep_cols], n_filtered[keep_cols]], ignore_index=True)
    return out


# -----------------------------------------------------------------------------
# CLI
# -----------------------------------------------------------------------------

@app.command("assign")
def assign_command(
    cluster_inputs: List[Path] = typer.Option(..., "-c", "--clusters",
                                              help="Per-sample cluster files (from clustering or consensusCluster)"),
    sample_names: str = typer.Option(..., "-n", "--sample-names",
                                     help='Sample names, space-separated'),
    gff_file: Path = typer.Option(..., "-a", "--annotation",
                                  help="GFF3 annotation file"),
    output_prefix: Path = typer.Option(..., "-o", "--output-prefix",
                                       help="Output prefix; produces <prefix>.<sample>.{assigned,unassigned,filtered}.tsv"),
    upstream: int = typer.Option(1000, "--upstream"),
    upstream_overlap: int = typer.Option(500, "--upstream-overlap"),
    downstream: int = typer.Option(0, "--downstream"),
    filter_cluster: bool = typer.Option(True, "--filter-cluster/--no-filter-cluster"),
    filter_cluster_threshold: float = typer.Option(0.02, "--filter-cluster-threshold"),
    ref_table: Optional[Path] = typer.Option(None, "--ref-table",
                                             help="Optional precomputed gene-reference TSV (parity testing)"),
):
    """TSSr-style annotateCluster: assign clusters to genes via GFF."""
    names = sample_names.split()
    if len(names) != len(cluster_inputs):
        raise typer.BadParameter(
            f"sample-names count ({len(names)}) != cluster files count ({len(cluster_inputs)})")

    if ref_table is not None:
        ref_genes = pd.read_csv(ref_table, sep='\t')
    else:
        ref_genes = load_genes_from_gff(str(gff_file))
    logger.info(f"loaded {len(ref_genes)} gene features")

    for n, p in zip(names, cluster_inputs):
        cluster_df = pd.read_csv(p, sep='\t')
        asn = assign_one_sample(cluster_df, ref_genes,
                                upstream=upstream,
                                upstream_overlap=upstream_overlap,
                                downstream=downstream,
                                filter_cluster=filter_cluster)
        assigned = asn[asn['gene'].notna()].drop(columns=['inCoding'], errors='ignore').copy()
        unassigned = asn[asn['gene'].isna()].drop(columns=['inCoding'], errors='ignore').copy()
        assigned.to_csv(f"{output_prefix}.{n}.assigned.tsv", sep='\t', index=False)
        unassigned.to_csv(f"{output_prefix}.{n}.unassigned.tsv", sep='\t', index=False)
        logger.info(f"{n}: assigned={len(assigned)}  unassigned={len(unassigned)}")
        if filter_cluster:
            filtered = filter_assigned(asn, filter_cluster_threshold)
            filtered.to_csv(f"{output_prefix}.{n}.filtered.tsv", sep='\t', index=False)
            logger.info(f"{n}: filtered={len(filtered)}")


if __name__ == '__main__':
    app()
