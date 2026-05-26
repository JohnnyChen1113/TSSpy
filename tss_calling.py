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

__version__ = "0.13.0"  # --strict-tssr flag: true keeps v0.12.0 TSSr parity, false fixes two TSSr bugs

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


def remove_g_mismatch_with_reference(read, fasta: pysam.FastaFile, chrom: str,
                                     strict_tssr: bool = True) -> int:
    """
    Remove mismatched G at 5' end using reference genome sequence.

    Two modes via `strict_tssr`:
      True  — bit-for-bit reproduction of TSSr 0.99.6 getTSS(). On minus strand
              this uses TSSr's sum-all-CIGAR end (1bp past actual alignment for
              insertion reads) and the iteration filter checks seq[0:i] (the
              non-TSS end of the forward-orient seq). These are TSSr bugs that
              TSSpy reproduces faithfully in this mode.
      False — corrected behavior for biologically faithful TSS calls:
              minus-strand end uses pysam.reference_end (true alignment end),
              and the iteration filter checks seq[-i:] (the actual 5' end of
              the cDNA in BAM forward orientation).
    """
    read_seq = read.query_sequence
    if not read_seq:
        return None

    if read.is_reverse:
        read_len = len(read_seq)
        if strict_tssr:
            # TSSr-parity: end = BAM_POS + sum(all CIGAR ints) - 1
            pos = (read.reference_start + 1) + calculate_tssr_mapped_length(read.cigarstring) - 1
        else:
            # Corrected: true alignment end (reference-consuming CIGAR ops only)
            pos = read.reference_end  # 0-based exclusive == 1-based end

        # Round 1: last base of forward-orient seq must be C (= complement of cDNA 5' G)
        last_base = read_seq[read_len - 1]
        if last_base != 'C':
            return pos

        # Reference check at the current end position (1-based 'pos' → 0-based pos-1)
        ref_end_pos = pos - 1
        try:
            ref_base = fasta.fetch(chrom, ref_end_pos, ref_end_pos + 1).upper()
        except:
            return pos

        complement = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G', 'N': 'N'}
        ref_base_complement = complement.get(ref_base, 'N')

        if ref_base_complement == 'G':
            return pos

        pos -= 1

        # Iterative rounds: TSSr-parity checks seq[0:i] (wrong end, buggy);
        # corrected checks seq[-i:] (the actual cDNA 5' i bases in forward orientation).
        i = 2
        while i <= read_len:
            if strict_tssr:
                prefix = read_seq[0:i]
            else:
                prefix = read_seq[-i:]
            expected = 'C' * i
            if prefix != expected:
                break

            new_ref_pos = pos - 1
            try:
                ref_base = fasta.fetch(chrom, new_ref_pos, new_ref_pos + 1).upper()
            except:
                break

            ref_base_complement = complement.get(ref_base, 'N')
            if ref_base_complement == 'G':
                break

            pos -= 1
            i += 1

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


def get_tss_position_no_reference(read, strict_tssr: bool = True) -> tuple:
    """
    Get TSS position when no reference is provided.

    Returns:
        (position, strand) tuple
    """
    if read.is_reverse:
        if strict_tssr:
            pos = (read.reference_start + 1) + calculate_tssr_mapped_length(read.cigarstring) - 1
        else:
            pos = read.reference_end  # true alignment end
        strand = "-"
    else:
        pos = read.reference_start + 1  # Convert to 1-based
        strand = "+"

    return pos, strand


def process_single_bam(bam_file: str,
                       sample_name: str,
                       sequencing_quality_threshold: int,
                       mapping_quality_threshold: int,
                       reference_file: Optional[str] = None,
                       strict_tssr: bool = True) -> pd.DataFrame:
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

                # TSSr passes isNotPassingQualityControls=FALSE to scanBamFlag,
                # which excludes reads with the BAM 0x200 (QC fail) bit set.
                if read.is_qcfail:
                    continue

                if read.mapping_quality < mapping_quality_threshold:
                    continue

                if read.query_qualities:
                    avg_qual = sum(read.query_qualities) / len(read.query_qualities)
                    if avg_qual < sequencing_quality_threshold:
                        continue

                chrom = bam.get_reference_name(read.reference_id)

                if fasta:
                    pos = remove_g_mismatch_with_reference(read, fasta, chrom, strict_tssr=strict_tssr)
                    if pos is None:
                        continue
                    strand = "-" if read.is_reverse else "+"
                else:
                    pos, strand = get_tss_position_no_reference(read, strict_tssr=strict_tssr)

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
        10,
        "--sequencing-quality",
        help="Minimum sequencing quality threshold (TSSr default: 10)",
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
    strict_tssr: bool = typer.Option(
        True,
        "--strict-tssr/--no-strict-tssr",
        help=(
            "true (default): bit-for-bit reproduction of TSSr 0.99.6, including "
            "two known TSSr bugs (sum-all-CIGAR mapped_length; minus-strand iteration "
            "filter checking the wrong end of the read). "
            "false: apply both corrections — minus-strand end uses true alignment "
            "end and iteration checks the actual cDNA 5' bases. "
            "Diverges from TSSr output but is biologically faithful."
        ),
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
    print(f"Mode: {'strict TSSr parity' if strict_tssr else 'CORRECTED (--no-strict-tssr)'}")

    if processes is None:
        processes = max(1, cpu_count() - 1)
    processes = min(processes, len(input_files_str))
    print(f"Using {processes} processes")

    args_list = [
        (bam, name, sequencing_quality, mapping_quality, reference_str, strict_tssr)
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
