#!/usr/bin/env python3
"""
TSS Calling - Extract TSS information from BAM files
Supports reference-based G mismatch removal for CAGE data
"""

import typer
import pysam
import pandas as pd
from collections import defaultdict
from typing import List, Optional
import logging
from multiprocessing import Pool, cpu_count
from pathlib import Path
import re

# Configure logging - default to WARNING level
logging.basicConfig(level=logging.WARNING, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

__version__ = "0.11.0"  # Use reference_end for correct position (handles indels correctly)

app = typer.Typer(help=f"Extract TSS information from BAM files (v{__version__})")


def calculate_tssr_mapped_length(cigar_string: str) -> int:
    """
    Calculate mapped length the way TSSr does: sum of ALL CIGAR numbers.

    TSSr extracts all numbers from CIGAR and sums them, regardless of operation type.
    Example: "39M1I35M" -> 39 + 1 + 35 = 75
    Example: "32M1D42M" -> 32 + 1 + 42 = 75

    This differs from:
    - query_length: only counts operations that consume query (M, I, S, =, X)
    - reference_length: only counts operations that consume reference (M, D, N, =, X)
    """
    if not cigar_string:
        return 0
    numbers = re.findall(r'(\d+)', cigar_string)
    return sum(int(n) for n in numbers)


def remove_g_mismatch_with_reference(read, fasta: pysam.FastaFile, chrom: str) -> int:
    """
    Remove mismatched G at 5' end using reference genome sequence.

    Implements TSSr algorithm:
    - Plus strand: check if read starts with G and reference is not G
    - Minus strand: check if read ends with C (complement of G) and reference at START is not G

    In BAM format, minus strand reads are stored as reverse complement.
    RNA 5' cap G appears as C at the END of query_sequence.

    Args:
        read: pysam AlignedSegment
        fasta: pysam FastaFile for reference genome
        chrom: chromosome name

    Returns:
        Adjusted TSS position (1-based)
    """
    read_seq = read.query_sequence
    if not read_seq:
        return None

    if read.is_reverse:
        # Minus strand: TSS is at the "end" position (5' end of minus strand read)
        # Use pysam's reference_end which correctly handles CIGAR operations:
        # - M, D, N, =, X consume reference
        # - I, S, H, P do NOT consume reference
        # pysam's reference_end is 0-based exclusive, which equals the 1-based end position
        read_len = len(read_seq)
        pos = read.reference_end  # 0-based exclusive = 1-based end position

        # Minus strand G mismatch removal algorithm (following TSSr exactly):
        #
        # Key insight from R testing:
        # 1. BAM stores minus strand reads in FORWARD strand direction (NOT reverse complement!)
        # 2. TSSr uses resize(width=1, fix='start') which for minus strand returns END position
        # 3. TSSr uses getSeq() which returns COMPLEMENT of reference for minus strand
        # 4. TSSr checks if complement != 'G' (i.e., forward strand ref != 'C')
        #
        # Algorithm:
        # 1. Check if seq[-1] == 'C' (last base, corresponds to END position)
        # 2. Check reference at END position, get complement
        # 3. If complement != 'G', it's a mismatch -> end -= 1
        # 4. Iterate: check if seq[0:i] == 'CC...C' (first i bases are all C)

        # Round 1: Check if last base is C
        last_base = read_seq[read_len - 1]
        if last_base != 'C':
            return pos

        # Check reference at END position
        # pos is pysam's reference_end (0-based exclusive), so the last aligned base is at pos - 1
        ref_end_pos = pos - 1  # 0-based position of last aligned base
        try:
            ref_base = fasta.fetch(chrom, ref_end_pos, ref_end_pos + 1).upper()
        except:
            return pos

        # TSSr's getSeq() returns complement for minus strand
        # Check if complement != 'G' (i.e., forward strand ref != 'C')
        complement = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G', 'N': 'N'}
        ref_base_complement = complement.get(ref_base, 'N')

        if ref_base_complement == 'G':
            # Complement is G (i.e., ref is C), not a mismatch
            return pos

        # Mismatch! Decrement TSS
        pos -= 1

        # Iterative rounds: check for consecutive C at sequence START
        # TSSr checks: substr(seq, start=1, stop=i) == paste(rep("C",i), collapse="")
        # This checks if the FIRST i bases are all 'C'
        i = 2
        while i <= read_len:
            # Check if first i bases of sequence are all C
            prefix = read_seq[0:i]
            expected = 'C' * i
            if prefix != expected:
                break

            # Check reference at NEW end position (pos has been decremented)
            # The new end position is pos - 1 (0-based)
            new_ref_pos = pos - 1  # 0-based
            try:
                ref_base = fasta.fetch(chrom, new_ref_pos, new_ref_pos + 1).upper()
            except:
                break

            ref_base_complement = complement.get(ref_base, 'N')
            if ref_base_complement == 'G':
                break

            pos -= 1
            i += 1

            # Safety limit
            if i > 10:
                break

        return pos

    else:
        # Plus strand: TSS is at reference_start + 1 (1-based)
        pos = read.reference_start + 1

        # Check for consecutive G at read start
        removed_count = 0
        i = 0

        while i < len(read_seq):
            read_base = read_seq[i]
            if read_base != 'G':
                break

            # Check reference at current position
            ref_pos = read.reference_start + i  # 0-based
            try:
                ref_base = fasta.fetch(chrom, ref_pos, ref_pos + 1).upper()
            except:
                break

            if ref_base == 'G':
                # Reference is G, not a mismatch
                break

            removed_count += 1
            i += 1

            if removed_count > 10:
                break

        if removed_count > 3:
            logger.debug(f"More than 3 G removed at {chrom}:{pos} (strand: +)")

        # TSSr does: start += removed_count
        return pos + removed_count


def get_tss_position_no_reference(read) -> tuple:
    """
    Get TSS position when no reference is provided.

    Returns:
        (position, strand) tuple
    """
    if read.is_reverse:
        # Use reference_end (0-based exclusive = 1-based end position)
        pos = read.reference_end
        strand = "-"
    else:
        pos = read.reference_start + 1  # Convert to 1-based
        strand = "+"

    return pos, strand


def process_single_bam(bam_file: str,
                       sample_name: str,
                       sequencing_quality_threshold: int,
                       mapping_quality_threshold: int,
                       reference_file: Optional[str] = None) -> pd.DataFrame:
    """
    Process a single BAM file to extract TSS information.
    """
    try:
        fasta = None
        if reference_file:
            try:
                fasta = pysam.FastaFile(reference_file)
            except Exception as e:
                pass

        with pysam.AlignmentFile(bam_file, "rb") as bam:
            tss_dict = defaultdict(int)

            for read in bam:
                if read.is_unmapped:
                    continue

                if read.mapping_quality < mapping_quality_threshold:
                    continue

                if read.query_qualities:
                    avg_qual = sum(read.query_qualities) / len(read.query_qualities)
                    if avg_qual < sequencing_quality_threshold:
                        continue

                chrom = bam.get_reference_name(read.reference_id)

                if fasta:
                    pos = remove_g_mismatch_with_reference(read, fasta, chrom)
                    if pos is None:
                        continue
                    strand = "-" if read.is_reverse else "+"
                else:
                    pos, strand = get_tss_position_no_reference(read)

                tss_dict[(chrom, pos, strand)] += 1

            if fasta:
                fasta.close()

            if tss_dict:
                df = pd.DataFrame(
                    [(k[0], k[1], k[2], v) for k, v in tss_dict.items()],
                    columns=['chr', 'pos', 'strand', 'count']
                )
                df = df.rename(columns={'count': sample_name})
                return df

            return pd.DataFrame(columns=['chr', 'pos', 'strand'])

    except Exception as e:
        import traceback
        traceback.print_exc()
        return pd.DataFrame(columns=['chr', 'pos', 'strand'])


def _process_bam_wrapper(args):
    """Wrapper function for multiprocessing"""
    return process_single_bam(*args)


@app.callback(invoke_without_command=True)
def main(
    ctx: typer.Context,
    input_files: List[Path] = typer.Option(
        None,
        "-i", "--input",
        exists=True,
        help="Input BAM files (can specify multiple times)",
    ),
    output_file: Path = typer.Option(
        None,
        "-o", "--output",
        help="Output TSS table file",
    ),
    sample_names: Optional[str] = typer.Option(
        None,
        "-n", "--sample-names",
        help="Sample names (space-separated)",
    ),
    reference: Optional[Path] = typer.Option(
        None,
        "-r", "--reference",
        help="Reference genome FASTA file (for G mismatch removal)",
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
    processes: Optional[int] = typer.Option(
        None,
        "-p", "--processes",
        help="Number of processes (default: CPU cores - 1)",
    ),
    verbose: bool = typer.Option(
        False,
        "-v", "--verbose",
        help="Enable verbose output",
    ),
):
    """
    Extract TSS information from BAM files.

    Example:
        tsspy tssCalling -i S01.bam -i S02.bam -o output.tsv -n "sample1 sample2" -r genome.fa
    """
    if ctx.invoked_subcommand is not None:
        return

    if not input_files or not output_file:
        print(ctx.get_help())
        raise typer.Exit(0)

    if verbose:
        logging.getLogger().setLevel(logging.INFO)

    input_files_str = [str(f) for f in input_files]
    output_file_str = str(output_file)
    reference_str = str(reference) if reference else None

    if sample_names is not None:
        names_list = sample_names.split()
        if len(names_list) != len(input_files_str):
            raise typer.BadParameter(
                f"Number of sample names ({len(names_list)}) must match input files ({len(input_files_str)})"
            )
    else:
        names_list = [Path(f).stem for f in input_files_str]

    print(f"Processing {len(input_files_str)} BAM files")
    print(f"Sample names: {', '.join(names_list)}")
    if reference_str:
        print(f"G mismatch removal: enabled")
    else:
        print(f"G mismatch removal: disabled")

    if processes is None:
        processes = max(1, cpu_count() - 1)
    processes = min(processes, len(input_files_str))
    print(f"Using {processes} processes")

    args_list = [
        (bam, name, sequencing_quality, mapping_quality, reference_str)
        for bam, name in zip(input_files_str, names_list)
    ]

    if processes > 1:
        with Pool(processes=processes) as pool:
            results = pool.map(_process_bam_wrapper, args_list)
    else:
        results = [_process_bam_wrapper(args) for args in args_list]

    all_samples = [df for df in results if not df.empty]

    if not all_samples:
        print("Error: No valid data found")
        raise typer.Exit(1)

    print("Merging samples...")
    final_df = all_samples[0]
    for df in all_samples[1:]:
        final_df = pd.merge(final_df, df, on=['chr', 'pos', 'strand'], how='outer')

    final_df = final_df.fillna(0)
    for col in final_df.columns:
        if col not in ['chr', 'pos', 'strand']:
            final_df[col] = final_df[col].astype(int)

    print("Sorting...")
    strand_order = {'+': 0, '-': 1}
    final_df['strand_order'] = final_df['strand'].map(strand_order)
    final_df = final_df.sort_values(['strand_order', 'chr', 'pos'])
    final_df = final_df.drop('strand_order', axis=1)

    print(f"Saving to {output_file_str}")
    final_df.to_csv(output_file_str, sep='\t', index=False)
    print(f"Done! Total: {len(final_df)} positions")


if __name__ == '__main__':
    app()
