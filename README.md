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
pip install typer[all] pysam pandas numpy biopython pyBigWig
```

### Using conda (recommended)

```bash
conda create -n tsspy python=3.8
conda activate tsspy
conda install -c bioconda -c conda-forge pysam pandas numpy typer biopython pybigwig
```

## Dependencies
- Python >= 3.8
- typer
- pysam
- pandas
- numpy
- biopython
- pyBigWig

## Quick Start

### 1. TSS Calling
Extract TSS from BAM files:
```bash
python -m TSSpy.main tssCalling -i sample1.bam -o sample1.TSS.tsv
```
Multiple BAM files in parallel:
```bash
python -m TSSpy.main tssCalling -i sample1.bam -i sample2.bam -o all.TSS.tsv -n "sample1 sample2"
```

### 2. TSS Clustering
Cluster TSSs to infer core promoters:
```bash
python -m TSSpy.main clustering -i all.TSS.tsv -o all.TSS.clustered.tsv -s sample1
```

### 3. Generate BigWig/BedGraph Files
Create visualization files for genome browsers:
```bash
# Generate bigWig files (requires chromosome sizes)
python -m TSSpy.main bigwig --input all.TSS.tsv --output-prefix all_samples --format bigwig --reference genome.fa --process

# Generate bedGraph files
python -m TSSpy.main bigwig --input all.TSS.tsv --output-prefix all_samples --format bedgraph
```

### 4. Correlation Analysis
Calculate correlations between samples:
```bash
python -m TSSpy.main correlation --input all.TSS.tsv --output correlation_matrix.tsv
```

### 5. Plot TSS Data
Generate visualization plots:
```bash
python -m TSSpy.main plot --help
```

### 6. Gene Assignment (Under Development)
```bash
python -m TSSpy.main geneAssign --help
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
- Built with modern CLI frameworks (Typer/Click)

## Publishing to Bioconda

To publish TSSpy to Bioconda, follow these steps:

1. **Prepare your package:**
   - Ensure you have a proper `setup.py` file
   - Create/update your PyPI package
   - Make sure your package has a version number in `__init__.py`

2. **Fork the bioconda-recipes repository:**
   - Visit https://github.com/bioconda/bioconda-recipes
   - Create a fork to your GitHub account

3. **Create a new recipe:**
   - Clone your forked repository
   - Create a new directory in `recipes/tsspy/`
   - Add a `meta.yaml` file with package information:

```yaml
{% set name = "tsspy" %}
{% set version = "0.1.0" %}  # Update with your version

package:
  name: {{ name|lower }}
  version: {{ version }}

source:
  url: https://github.com/yourusername/TSSpy/archive/v{{ version }}.tar.gz
  sha256: # generate this with `sha256sum` or similar tool

build:
  number: 0
  noarch: python
  script: {{ PYTHON }} -m pip install . --no-deps --ignore-installed -vv

requirements:
  host:
    - python >=3.7
    - pip
    - setuptools
  run:
    - python >=3.7
    - typer
    - pysam
    - pandas
    - numpy
    - biopython
    - pybigwig

test:
  imports:
    - TSSpy
  commands:
    - python -m TSSpy.main --help

about:
  home: https://github.com/yourusername/TSSpy
  license: MIT
  summary: 'Python command-line tool for transcription start site (TSS) analysis'
  description: |
    TSSpy is a Python command-line tool designed for the analysis of 
    transcription start site (TSS) data. Inspired by the TSSr (R/Bioconductor) 
    package, TSSpy adopts a main-command + subcommand architecture for easy 
    extension and maintenance.
  dev_url: https://github.com/yourusername/TSSpy

extra:
  recipe-maintainers:
    - yourusername
```

4. **Test your recipe locally:**
   - Install conda-build: `conda install conda-build`
   - Build locally: `conda build recipes/tsspy`

5. **Submit a pull request:**
   - Commit your changes and push to your fork
   - Create a PR to the bioconda-recipes repository
   - Wait for the automated tests and review process

6. **After PR is merged:**
   - Your package will be built and available in the bioconda channel
   - Users can install with: `conda install -c bioconda tsspy`

For more detailed instructions, see the [Bioconda documentation](https://bioconda.github.io/contributor/index.html).

## Contributing
Contributions, issues, and feature requests are welcome!

## License
MIT 