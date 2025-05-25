import typer
import pandas as pd
import numpy as np
import pyBigWig
import pysam
from Bio import SeqIO
import tempfile
import os
from pathlib import Path
from typing import Optional

app = typer.Typer(help="Generate BigWig/BedGraph files from TSS table or BAM file.")

def extract_chrom_sizes_from_fasta(fasta_path: str, out_path: str):
    with open(out_path, 'w') as f:
        for record in SeqIO.parse(fasta_path, "fasta"):
            f.write(f"{record.id}\t{len(record.seq)}\n")

def extract_chrom_sizes_from_bam(bam_path: str, out_path: str):
    with pysam.AlignmentFile(bam_path, "rb") as bam:
        with open(out_path, 'w') as f:
            for ref, length in zip(bam.references, bam.lengths):
                f.write(f"{ref}\t{length}\n")

def filter_tss_table(df, process, abs_threshold, tpm_threshold):
    df_filtered = df.copy()
    sample_cols = df.columns[3:]
    if process:
        # Normalize to TPM
        for col in sample_cols:
            total = df_filtered[col].sum()
            if total > 0:
                df_filtered[col] = df_filtered[col] / total * 1e6
            else:
                df_filtered[col] = 0
        # Filter by TPM threshold
        mask = (df_filtered[sample_cols] >= tpm_threshold).any(axis=1)
        df_filtered = df_filtered[mask]
    else:
        # Filter by absolute threshold
        mask = (df_filtered[sample_cols] >= abs_threshold).any(axis=1)
        df_filtered = df_filtered[mask]
    return df_filtered

def write_bedgraph(df, sample, strand, out_path):
    # bedGraph: chr, start, end, value
    sub = df[df['strand'] == strand][['chr', 'pos', sample]]
    sub = sub[sub[sample] > 0]
    
    # Sort by chromosome and position
    sub = sub.sort_values(['chr', 'pos'])
    
    with open(out_path, 'w') as f:
        for _, row in sub.iterrows():
            chrom = row['chr']
            start = int(row['pos']) - 1
            end = int(row['pos'])
            value = row[sample]
            f.write(f"{chrom}\t{start}\t{end}\t{value}\n")

def write_bigwig(df, sample, strand, chrom_sizes_path, out_path):
    # Load chromosome sizes
    chrom_sizes = []
    with open(chrom_sizes_path) as f:
        for line in f:
            chrom, size = line.strip().split('\t')
            chrom_sizes.append((chrom, int(size)))
    
    # Get chromosome order
    chr_order = {chrom: i for i, (chrom, _) in enumerate(chrom_sizes)}
    
    # Filter data and add sorting column
    sub = df[df['strand'] == strand][['chr', 'pos', sample]].copy()
    sub = sub[sub[sample] > 0]
    
    # Ensure chromosomes exist in chrom_sizes
    valid_chroms = set(chr_order.keys())
    sub = sub[sub['chr'].isin(valid_chroms)]
    
    if sub.empty:
        typer.echo(f"Warning: No data for {sample} on {strand} strand after filtering.")
        # Create empty bigWig file
        bw = pyBigWig.open(out_path, "w")
        bw.addHeader(chrom_sizes)
        bw.close()
        return
    
    # Add chromosome order field
    sub['chr_order'] = sub['chr'].map(chr_order)
    
    # Sort by chromosome order and position
    sub = sub.sort_values(['chr_order', 'pos'])
    
    # Initialize bigWig file
    bw = pyBigWig.open(out_path, "w")
    bw.addHeader(chrom_sizes)
    
    # Process by chromosome to ensure proper ordering
    for chrom, _ in chrom_sizes:
        chrom_data = sub[sub['chr'] == chrom]
        if chrom_data.empty:
            continue
            
        # Sort by position
        chrom_data = chrom_data.sort_values('pos')
        
        # Prepare data for this chromosome
        chroms = [chrom] * len(chrom_data)
        starts = (chrom_data['pos'] - 1).astype(int).tolist()
        ends = chrom_data['pos'].astype(int).tolist()
        values = chrom_data[sample].astype(float).tolist()
        
        # Ensure all lists have the same length
        if len(starts) != len(ends) or len(starts) != len(values):
            typer.echo(f"Warning: Length mismatch for {chrom}: starts={len(starts)}, ends={len(ends)}, values={len(values)}")
            continue
        
        # Add entries for this chromosome
        if starts:
            try:
                bw.addEntries(chroms, starts, ends=ends, values=values)
            except Exception as e:
                typer.echo(f"Error adding entries for {chrom}: {e}")
    
    bw.close()

@app.callback(invoke_without_command=True)
def main(ctx: typer.Context,
    input: Path = typer.Option(..., "-i", "--input", help="Input TSS table (tsv) or BAM file"),
    output_prefix: str = typer.Option(..., "-o", "--output-prefix", help="Output file prefix"),
    format: str = typer.Option("bigwig", "--format", help="Output format: bigwig or bedgraph"),
    process: bool = typer.Option(False, "--process", help="Normalize to TPM before exporting"),
    raw: bool = typer.Option(False, "--raw", help="Use raw counts (default)"),
    abs_threshold: float = typer.Option(2, "--abs-threshold", help="Absolute value threshold (default: 2)"),
    tpm_threshold: float = typer.Option(0.02, "--tpm-threshold", help="TPM threshold (default: 0.02)"),
    reference: Optional[Path] = typer.Option(None, "--reference", help="Reference genome fasta (for chrom.sizes)"),
    bam: Optional[Path] = typer.Option(None, "--bam", help="BAM file (for chrom.sizes)"),
    chrom_sizes: Optional[Path] = typer.Option(None, "--chrom-sizes", help="Existing chrom.sizes file")
):
    """
    Generate BigWig/BedGraph files from TSS table or BAM file.
    """
    # Skip execution if called as a group
    if ctx.invoked_subcommand is not None:
        return
        
    # 1. Read TSS table or BAM
    if str(input).endswith('.bam'):
        typer.echo("[ERROR] Direct conversion from BAM to bigwig not implemented yet. Please use tssCalling first!")
        raise typer.Exit(1)
    else:
        try:
            df = pd.read_csv(input, sep='\t')
            typer.echo(f"Loaded TSS table with {len(df)} rows and {len(df.columns)} columns.")
        except Exception as e:
            typer.echo(f"[ERROR] Failed to read TSS table: {e}")
            raise typer.Exit(1)
    
    # Check table format
    required_cols = ['chr', 'pos', 'strand']
    if not all(col in df.columns for col in required_cols):
        typer.echo(f"[ERROR] TSS table must contain columns: {', '.join(required_cols)}")
        raise typer.Exit(1)
    
    # 2. Normalize and filter
    df = filter_tss_table(df, process, abs_threshold, tpm_threshold)
    sample_cols = [col for col in df.columns if col not in ['chr', 'pos', 'strand']]
    
    if not sample_cols:
        typer.echo("[ERROR] No sample columns found in the table!")
        raise typer.Exit(1)
        
    typer.echo(f"Found {len(sample_cols)} sample columns: {', '.join(sample_cols)}")
    
    # 3. Handle chrom.sizes
    need_chrom_sizes = (format.lower() == 'bigwig')
    chrom_sizes_path = None
    tmp_chrom_sizes = None
    
    if need_chrom_sizes:
        if chrom_sizes:
            chrom_sizes_path = chrom_sizes
            typer.echo(f"Using provided chromosome sizes file: {chrom_sizes}")
        elif reference:
            tmp_chrom_sizes = tempfile.NamedTemporaryFile(delete=False, suffix='.chrom.sizes')
            extract_chrom_sizes_from_fasta(str(reference), tmp_chrom_sizes.name)
            chrom_sizes_path = tmp_chrom_sizes.name
            typer.echo(f"Extracted chromosome sizes from reference genome: {reference}")
        elif bam:
            tmp_chrom_sizes = tempfile.NamedTemporaryFile(delete=False, suffix='.chrom.sizes')
            extract_chrom_sizes_from_bam(str(bam), tmp_chrom_sizes.name)
            chrom_sizes_path = tmp_chrom_sizes.name
            typer.echo(f"Extracted chromosome sizes from BAM file: {bam}")
        else:
            typer.echo("[ERROR] For bigwig output, you must specify --reference, --bam, or --chrom-sizes!")
            raise typer.Exit(1)
    
    # 4. Export files
    for sample in sample_cols:
        for strand, strand_label in [('+', 'plus'), ('-', 'minus')]:
            out_path = f"{output_prefix}.{sample}.{strand_label}.{format.lower() if format.lower() == 'bedgraph' else 'bw'}"
            
            if format.lower() == 'bedgraph':
                write_bedgraph(df, sample, strand, out_path)
                typer.echo(f"[bedGraph] {out_path}")
            elif format.lower() == 'bigwig':
                try:
                    write_bigwig(df, sample, strand, chrom_sizes_path, out_path)
                    typer.echo(f"[bigWig] {out_path}")
                except Exception as e:
                    typer.echo(f"[ERROR] Failed to create bigWig file {out_path}: {e}")
                    continue
    
    # 5. Clean up temporary files
    if tmp_chrom_sizes:
        os.unlink(tmp_chrom_sizes.name)
        typer.echo("Temporary chromosome sizes file removed.")

if __name__ == "__main__":
    app() 