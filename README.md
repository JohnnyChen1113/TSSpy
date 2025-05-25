# TSSpy: Python CLI for TSS Analysis

TSSpy is a Python command-line tool designed for the analysis of transcription start site (TSS) data. Inspired by the TSSr (R/Bioconductor) package (https://github.com/Linlab-slu/TSSr)

## Features
- Main-command + subcommand structure for modularity and scalability
- Multi-processing support for fast analysis
- Compatible with BAM and TSS table input formats
- Supports BigWig/BedGraph export for genome browser visualization
- Clustering, correlation analysis, and promoter identification
- Suitable for large-scale high-throughput sequencing data
- Clean codebase for easy secondary development

## Installation

### Using pip
It is recommended to use Python 3.8 or above.

```bash
pip install tsspy
```

After installation, you can use the `tsspy` command directly from the terminal.

### Using conda (recommended)

```bash
conda create -n tsspy python=3.8
conda activate tsspy
conda install -c bioconda -c conda-forge tsspy
```

## Dependencies
- Python >= 3.8
- typer (automatically installs click as a dependency)
- pysam
- pandas
- numpy
- biopython
- pyBigWig

## Quick Start

### Command Usage

After installation, you can use TSSpy in two ways:

1. **Direct command (recommended):**
   ```bash
   tsspy [command] [options]
   ```

2. **Module invocation (for development):**
   ```bash
   python -m TSSpy.main [command] [options]
   ```

All examples below use the direct command format.

### 1. TSS Calling
Extract TSS from BAM files:
```bash
tsspy tssCalling -i sample1.bam -o sample1.TSS.tsv
```
Multiple BAM files in parallel:
```bash
tsspy tssCalling -i sample1.bam -i sample2.bam -o all.TSS.tsv -n "sample1 sample2"
```

### 2. TSS Clustering
Cluster TSSs to infer core promoters:
```bash
tsspy clustering -i all.TSS.tsv -o all.TSS.clustered.tsv -s sample1
```

### 3. Generate BigWig/BedGraph Files
Create visualization files for genome browsers:
```bash
# Generate bigWig files (requires chromosome sizes)
tsspy bigwig --input all.TSS.tsv --output-prefix all_samples --format bigwig --reference genome.fa --process

# Generate bedGraph files
tsspy bigwig --input all.TSS.tsv --output-prefix all_samples --format bedgraph
```

### 4. Correlation Analysis
Calculate correlations between samples:
```bash
tsspy correlation --input all.TSS.tsv --output correlation_matrix.tsv
```

### 5. Plot TSS Data
Generate visualization plots:
```bash
tsspy plot --help
```

### 6. Gene Assignment (Under Development)
```bash
tsspy geneAssign --help
```

## Subcommands
- `tssCalling`: Identify TSSs from BAM files and output a TSS table
- `clustering`: Cluster TSSs to infer core promoters
- `bigwig`: Generate BigWig/BedGraph files for genome browser visualization
- `correlation`: Calculate correlations between TSS samples
- `plot`: Generate plots and visualizations for TSS data
- `geneAssign`: Assign TSS clusters to genes (under development)

## Extensibility
- To add a new feature, simply create a new subcommand module in `TSSpy/` and register it in `main.py`
- Shared utility functions can be placed in `TSSpy/utils.py`
- Built with modern CLI frameworks

## Contributing
Contributions, issues, and feature requests are welcome!

## License
MIT 