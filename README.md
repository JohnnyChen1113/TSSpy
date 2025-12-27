# TSSpy: Python CLI for TSS Analysis

TSSpy is a Python command-line tool designed for the analysis of transcription start site (TSS) data. Inspired by [TSSr](https://github.com/Linlab-slu/TSSr) (R/Bioconductor) but implemented in pure Python **without BSgenome dependency**.

## Features

- **TSS Calling**: Extract TSS from BAM files with reference-based G mismatch removal
- **Sample Merging**: Merge biological replicates and normalize to TPM
- **TSS Clustering**: Cluster TSS to infer core promoters (peakclu algorithm)
- **Consensus Clustering**: Cross-sample consensus cluster aggregation
- **Promoter Shape Analysis**: Calculate PSS (Promoter Shape Score) and SI (Shape Index)
- **Gene Assignment**: Assign clusters to genes using GTF/GFF annotations
- **Visualization**: BigWig/BedGraph export, correlation plots
- **Multi-processing**: Parallel processing support for large datasets
- **Clean CLI**: Main-command + subcommand structure for modularity

## Installation

### Using pip (recommended)

```bash
pip install tsspy
```

### Using conda

```bash
conda create -n tsspy python=3.8
conda activate tsspy
conda install -c bioconda -c conda-forge tsspy
```

### From source

```bash
git clone https://github.com/JohnnyChen1113/TSSpy.git
cd TSSpy/TSSpy
pip install -e .
```

## Dependencies

- Python >= 3.8
- typer
- pysam
- pandas
- numpy
- biopython
- pyBigWig
- matplotlib
- scipy (for correlation plots)

## Quick Start

### Command Structure

```bash
tsspy <command> [subcommand] [options]
```

### Available Commands

| Command | Description |
|---------|-------------|
| `tssCalling` | Extract TSS from BAM files |
| `mergeSamples` | Merge samples and normalize data |
| `clustering` | Cluster TSS to infer core promoters |
| `consensusCluster` | Create consensus clusters across samples |
| `shapeCluster` | Calculate promoter shape scores (PSS/SI) |
| `geneAssign` | Assign clusters to genes |
| `bigwig` | Generate BigWig/BedGraph files |
| `correlation` | Calculate sample correlations |
| `plot` | Generate visualization plots |

## Workflow Example

### 1. TSS Calling from BAM files

```bash
# Basic usage (soft-clipping mode)
tsspy tssCalling main -i S01.bam -i S02.bam -o raw.TSS.tsv -n "sample1 sample2"

# With reference genome (G mismatch removal enabled)
tsspy tssCalling main -i S01.bam -i S02.bam -o raw.TSS.tsv \
    -n "sample1 sample2" -r reference.fa
```

### 2. Merge Samples and Normalize

```bash
# Merge biological replicates
tsspy mergeSamples merge -i raw.TSS.tsv -o merged.TSS.tsv \
    -g "control treat" -m "1 1 2 2"

# Normalize to TPM
tsspy mergeSamples normalize -i merged.TSS.tsv -o normalized.TSS.tsv

# One-step processing (merge + normalize + filter)
tsspy mergeSamples process -i raw.TSS.tsv -o processed.TSS.tsv \
    -g "control treat" -m "1 1 2 2" \
    --normalize --filter tpm --filter-threshold 0.1
```

### 3. TSS Clustering

```bash
tsspy clustering main -i processed.TSS.tsv -o clusters.tsv -s control \
    --peak-distance 100 --extension-distance 30 --cluster-threshold 1
```

### 4. Consensus Clustering

```bash
# From per-sample cluster files
tsspy consensusCluster cluster \
    -i control.clusters.tsv -i treat.clusters.tsv \
    -n "control treat" -o consensus.tsv -d 50

# Directly from TSS table
tsspy consensusCluster from-tss -t processed.TSS.tsv -o consensus.tsv
```

### 5. Promoter Shape Analysis (PSS)

```bash
# Calculate PSS for a sample
tsspy shapeCluster calculate -c clusters.tsv -t processed.TSS.tsv \
    -o shape.tsv -s control -m PSS

# Calculate for all samples
tsspy shapeCluster batch -c clusters.tsv -t processed.TSS.tsv \
    -o shape -m PSS

# Classify promoters as sharp/broad
tsspy shapeCluster classify -i shape.tsv -o classified.tsv -m PSS
```

### 6. Gene Assignment

```bash
# Assign clusters to genes
tsspy geneAssign assign -c clusters.tsv -a annotation.gtf \
    -o assigned.tsv --upstream 1000

# Create gene-level summary
tsspy geneAssign summary -i assigned.tsv -o gene_summary.tsv

# Filter to primary promoters only
tsspy geneAssign filter -i assigned.tsv -o primary.tsv --keep-primary
```

### 7. Generate BigWig Files

```bash
tsspy bigwig -i processed.TSS.tsv -o output_prefix \
    --format bigwig --reference genome.fa --process
```

### 8. Correlation Analysis

```bash
tsspy correlation -i processed.TSS.tsv -o correlation.csv --plot
```

## Shape Score Methods

### PSS (Promoter Shape Score)
- **Formula**: `PSS = -sum(p_i * log2(p_i)) * log2(IQW)`
- Lower PSS = sharper promoter
- PSS = 0 for singletons
- Reference: Lu and Lin 2019

### SI (Shape Index)
- **Formula**: `SI = 2 + sum(p_i * log2(p_i))`
- Higher SI = sharper promoter
- SI = 2 for singletons
- Reference: Hoskins et al. 2011

## G Mismatch Removal

When using CAGE data, the 5' end of reads may contain non-coded G bases (m7G cap). TSSpy can remove these mismatched G bases using a reference genome:

```bash
tsspy tssCalling main -i sample.bam -o output.tsv -r reference.fa
```

This requires:
1. Reference FASTA file (with index `.fai`)
2. `--allow-softclipping` must NOT be set

## Output Formats

### TSS Table
```
chr    pos    strand    sample1    sample2    ...
chrI   1000   +         10         15
chrI   1050   +         5          8
```

### Cluster Table
```
cluster  chr   start  end   strand  dominant_tss  tags  tags_TPM  q_0.1  q_0.9  interquantile_width
1        chrI  1000   1100  +       1050          100   50.5      1010   1080   71
```

## Differences from TSSr

| Feature | TSSr (R) | TSSpy (Python) |
|---------|----------|----------------|
| Reference genome | BSgenome package required | FASTA file (pysam) |
| BAM processing | Rsamtools | pysam |
| BigWig export | rtracklayer | pyBigWig |
| CLI framework | R functions | typer (modern CLI) |
| Parallelization | parallel (R) | multiprocessing (Python) |

## Contributing

Contributions, issues, and feature requests are welcome!

## Citation

If you use TSSpy, please cite:
- TSSpy: https://github.com/JohnnyChen1113/TSSpy
- TSSr: Lu, Z., Berry, K., Hu, Z., Zhan, Y., Ahn, T., & Lin, Z. (2021). TSSr: an R package for comprehensive analyses of TSS sequencing data. NAR Genomics and Bioinformatics, 3(4).

## License

MIT
