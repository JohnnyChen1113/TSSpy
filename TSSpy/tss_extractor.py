#!/usr/bin/env python3
"""
TSS Extractor - Extract TSS information from BAM files
Version: 0.7.5
"""

import typer
import pysam
import pandas as pd
from collections import defaultdict
from typing import List, Dict, Tuple, Optional
import logging
from multiprocessing import Pool, cpu_count
from functools import partial
import os
import re
from pathlib import Path

# Configure logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

# Version information
__version__ = "0.2.0"

# Create Typer app
app = typer.Typer(help=f"Extract TSS information from BAM files and output a table (v{__version__})")

def get_mapped_length(cigar: str, softclipping_allowed: bool) -> int:
    """
    Calculate mapped length from CIGAR string, similar to R code's logic
    """
    if not cigar:
        return 0
        
    # Extract all numbers from CIGAR string
    numbers = [int(n) for n in re.findall(r'(\d+)', cigar)]
    total_length = sum(numbers)
    
    if softclipping_allowed:
        # Extract soft clip length if present
        soft_clip_match = re.search(r'(\d+)S', cigar)
        if soft_clip_match:
            total_length -= int(soft_clip_match.group(1))
            
    return total_length

def parse_md_tag(md_tag: str) -> List[Tuple[int, str]]:
    """
    Parse MD tag to get mismatch information
    Example: MD:Z:0C74 -> [(0, 'C')] means first base is C in reference
    """
    if not md_tag:
        return []
    
    # Remove 'MD:Z:' prefix if present
    md = md_tag.split(':', 2)[-1]
    
    # Split by numbers and letters
    parts = re.findall(r'(\d+|[A-Z]|\^[A-Z]+)', md)
    
    result = []
    pos = 0
    for part in parts:
        if part.isdigit():
            pos += int(part)
        elif part.startswith('^'):
            # Deletion in reference
            pos += len(part) - 1
        else:
            # Mismatch
            result.append((pos, part))
            pos += 1
    
    return result

def adjust_pos_for_cage_mismatch(pos, mismatches, strand, chrom, sample_name, read_id=None, read_seq=None):
    """
    For CAGE data, only if the 5' end of the read is G, check for consecutive mismatches at the 5' end (both strands), regardless of reference base.
    If more than 3 consecutive mismatched G, print a warning including read id if available.
    """
    if not mismatches or not read_seq:
        return pos
    count = 0
    seq_len = len(read_seq)
    if strand == '+':
        if read_seq[0] != 'G':
            return pos
        for i, (mpos, _) in enumerate(mismatches):
            if mpos != i:
                break
            if read_seq[i] == 'G':
                count += 1
            else:
                break
        if count > 0:
            pos += count
            if count > 3:
                logger.warning(f"More than 3 consecutive mismatched G at 5' end on {chrom} (sample: {sample_name}, strand: +, read: {read_id}). Please check your data before further analysis.")
    elif strand == '-':
        if read_seq[-1] != 'G':
            return pos
        for i, (mpos, _) in enumerate(mismatches):
            if mpos != i:
                break
            if read_seq[-(i+1)] == 'G':
                count += 1
            else:
                break
        if count > 0:
            pos -= count
            if count > 3:
                logger.warning(f"More than 3 consecutive mismatched G at 5' end on {chrom} (sample: {sample_name}, strand: -, read: {read_id}). Please check your data before further analysis.")
    return pos

def process_single_bam(bam_file: str,
                      sample_name: str,
                      sequencing_quality_threshold: int,
                      mapping_quality_threshold: int,
                      softclipping_allowed: bool) -> pd.DataFrame:
    """
    Process a single BAM file to extract TSS information
    
    Args:
        bam_file: Path to BAM file
        sample_name: Name to use for this sample in output
        sequencing_quality_threshold: Minimum sequencing quality score
        mapping_quality_threshold: Minimum mapping quality score
        softclipping_allowed: Whether to allow soft clipping
        
    Returns:
        DataFrame containing TSS information
    """
    logger.info(f"Processing file: {bam_file}")
    try:
        # Check if BAM index exists
        bai_file = f"{bam_file}.bai"
        if softclipping_allowed and not os.path.exists(bai_file):
            logger.info(f"Creating BAM index for {bam_file}")
            pysam.index(bam_file)
            
        with pysam.AlignmentFile(bam_file, "rb") as bam:
            tss_dict = defaultdict(int)
            
            for read in bam:
                # Skip unmapped reads
                if read.is_unmapped:
                    continue
                    
                # Skip low quality reads
                if read.mapping_quality < mapping_quality_threshold:
                    continue
                    
                # Get average sequencing quality
                if read.query_qualities:
                    avg_qual = sum(read.query_qualities) / len(read.query_qualities)
                    if avg_qual < sequencing_quality_threshold:
                        continue
                
                # Get chromosome and position
                chrom = bam.get_reference_name(read.reference_id)
                
                # Calculate mapped length
                mapped_length = get_mapped_length(read.cigarstring, softclipping_allowed)
                
                if read.is_reverse:
                    # Minus strand
                    pos = read.reference_end
                    strand = "-"
                    if softclipping_allowed and read.cigarstring:
                        if read.cigarstring.endswith('S'):
                            soft_clip_len = int(re.search(r'(\d+)S$', read.cigarstring).group(1))
                            pos -= soft_clip_len
                    md_tag = read.get_tag('MD') if read.has_tag('MD') else None
                    if md_tag and read.query_sequence:
                        mismatches = parse_md_tag(md_tag)
                        pos = adjust_pos_for_cage_mismatch(pos, mismatches, strand, chrom, sample_name, read.query_name, read.query_sequence)
                else:
                    # Plus strand
                    pos = read.reference_start + 1  # Convert to 1-based
                    strand = "+"
                    if softclipping_allowed and read.cigarstring:
                        if read.cigarstring.startswith('S'):
                            soft_clip_len = int(re.search(r'^(\d+)S', read.cigarstring).group(1))
                            pos += soft_clip_len
                    md_tag = read.get_tag('MD') if read.has_tag('MD') else None
                    if md_tag and read.query_sequence:
                        mismatches = parse_md_tag(md_tag)
                        pos = adjust_pos_for_cage_mismatch(pos, mismatches, strand, chrom, sample_name, read.query_name, read.query_sequence)
                
                # Add to dictionary
                tss_dict[(chrom, pos, strand)] += 1
                
            # Convert to DataFrame
            if tss_dict:
                df = pd.DataFrame([(k[0], k[1], k[2], v) for k, v in tss_dict.items()],
                                columns=['chr', 'pos', 'strand', 'count'])
                # Rename count column to sample name
                df = df.rename(columns={'count': sample_name})
                return df
            return pd.DataFrame(columns=['chr', 'pos', 'strand'])
            
    except Exception as e:
        logger.error(f"Error processing {bam_file}: {str(e)}")
        return pd.DataFrame(columns=['chr', 'pos', 'strand'])

@app.command()
def main(
    input_files: List[Path] = typer.Option(
        ...,
        "-i",
        "--input-files",
        exists=True,
        help="Input BAM files (can specify multiple times)",
    ),
    output_file: Path = typer.Option(
        ...,
        "-o",
        "--output-file",
        help="Output TSS table file",
    ),
    sample_names: str = typer.Option(
        None,
        "-n",
        "--sample-names",
        help="Sample names for output columns, space-separated (must match number of input files)",
    ),
    sequencing_quality: int = typer.Option(
        20,
        "--sequencing-quality",
        help="Minimum sequencing quality threshold",
    ),
    mapping_quality: int = typer.Option(
        20,
        "--mapping-quality",
        help="Minimum mapping quality threshold",
    ),
    allow_softclipping: bool = typer.Option(
        False,
        "--allow-softclipping/--no-softclipping",
        help="Whether to allow soft clipping (default: False)",
    ),
    processes: Optional[int] = typer.Option(
        None,
        "--processes",
        help="Number of processes to use (default: number of CPU cores)",
    ),
):
    """
    Extract TSS information from BAM files and output a table.
    Multiple input files can be specified using -i multiple times.
    Example: python tss_extractor.py -i file1.bam -i file2.bam -i file3.bam -o output.tsv -n "sample1 sample2 sample3"
    """
    # Convert Path objects to strings
    input_files = [str(f) for f in input_files]
    output_file = str(output_file)
    
    # Parse sample names if provided
    if sample_names is not None:
        sample_names = sample_names.split()
        if len(sample_names) != len(input_files):
            raise typer.BadParameter(
                f"Number of sample names ({len(sample_names)}) must match number of input files ({len(input_files)})"
            )
    else:
        # Use file names as sample names if not provided
        sample_names = [Path(f).stem for f in input_files]
    
    logger.info(f"Processing {len(input_files)} input files: {', '.join(input_files)}")
    logger.info(f"Using sample names: {', '.join(sample_names)}")
    
    # Determine number of processes to use
    if processes is None:
        processes = max(1, cpu_count() - 1)  # Leave one CPU core free
    processes = min(processes, len(input_files))  # Don't use more processes than files
    
    logger.info(f"Using {processes} processes for parallel processing")
    
    # Create a partial function with fixed parameters
    process_func = partial(
        process_single_bam,
        sequencing_quality_threshold=sequencing_quality,
        mapping_quality_threshold=mapping_quality,
        softclipping_allowed=allow_softclipping,
    )
    
    # Process files in parallel
    with Pool(processes=processes) as pool:
        results = pool.starmap(
            process_func,
            zip(input_files, sample_names)
        )
    
    # Filter out empty DataFrames
    all_samples = [df for df in results if not df.empty]
    
    if not all_samples:
        logger.error("No valid data found in any input files")
        raise typer.Exit(1)
        
    # Merge all samples
    logger.info("Merging all samples...")
    final_df = all_samples[0]
    for df in all_samples[1:]:
        final_df = pd.merge(final_df, df, on=['chr', 'pos', 'strand'], how='outer')
    
    # Fill NA with 0 and convert to integer
    final_df = final_df.fillna(0)
    for col in final_df.columns:
        if col not in ['chr', 'pos', 'strand']:
            final_df[col] = final_df[col].astype(int)
    
    # Sort by strand (plus first), then chromosome and position
    logger.info("Sorting results...")
    strand_order = {'+': 0, '-': 1}
    final_df['strand_order'] = final_df['strand'].map(strand_order)
    final_df = final_df.sort_values(['strand_order', 'chr', 'pos'])
    final_df = final_df.drop('strand_order', axis=1)
    
    # Save to file
    logger.info(f"Saving results to {output_file}")
    final_df.to_csv(output_file, sep='\t', index=False)
    logger.info("Done!")

if __name__ == '__main__':
    app() 