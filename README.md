# TSSpy

**A Python CLI for transcription-start-site (TSS) analysis — bit-for-bit port of [TSSr](https://github.com/Linlab-slu/TSSr) 0.99.6, without the BSgenome dependency.**

TSSpy reproduces TSSr's full pipeline (`getTSS` → `mergeSamples` → `normalizeTSS` → `filterTSS` → `clusterTSS` → `consensusCluster` → `shapeCluster` → `annotateCluster` → `exportTSStoBedgraph`) in Python with `pysam` instead of Rsamtools/BSgenome.

## Parity guarantee

Every stage is regression-tested against TSSr 0.99.6 output on the 4-BAM yeast validation set:

| Stage | Parity |
|---|---|
| `tssCalling` | bit-for-bit (163,203 TSS positions, all counts) |
| `mergeSamples` / `normalize` / `filter` | bit-for-bit (TPM, poisson + TPM filters) |
| `clustering` (peakclu + peakcluMax) | bit-for-bit, including `peakcluMax` from the `adamZhang_TSSr` fork |
| `consensusCluster` | bit-for-bit |
| `shapeCluster` (PSS + SI) | within FP precision (max\|Δ\| < 1e-12, non-score columns exact) |
| `geneAssign` | bit-for-bit, including standalone GFF3 parser |
| `bigwig` (bedGraph + BigWig) | byte-for-byte (bedGraph) / interval-for-interval (BigWig) |

Run `pytest tests/` to verify (requires the BAM + GFF + TSSr-derived ground-truth files in the working dir; ~94s).

## Two correctness modes

The default `--strict-tssr` mode (= `true`) reproduces TSSr 0.99.6 exactly, including two known minus-strand bugs (see `TSSr_fix_proposal.md` for advisor write-up):

1. `mapped.length` overcounts insertion (I) operations → minus-strand TSS shifted 1bp past the actual alignment end for ~6,857 insertion-containing reads (S288C 4-BAM).
2. G-mismatch iteration filter checks `seq[1..i]` (the cDNA 3' end on minus strand) instead of the cDNA 5' bases — affects ~104K reads.

Pass `--no-strict-tssr` to apply the corrected algorithms. Validated YR motif uplift:

| Dataset | Read-weighted YR rate (minus strand) |
|---|---|
| S288C 4-BAM | 82.69% → 84.77% (+2.08 pp) |
| S. uvarum 2-BAM | 49.91% → 53.64% (+3.72 pp) |

## Install

```bash
git clone https://github.com/JohnnyChen1113/TSSpy.git
cd TSSpy
pip install -e .
```

This installs the `tsspy` console script and registers the `TSSpy` package.

Dependencies (auto-installed): `typer`, `pysam`, `pandas`, `numpy`, `biopython`, `pyBigWig`, `matplotlib`, `scipy`.

## Quick start

```bash
# 1. Call TSSs from BAMs (TSSr-strict by default)
tsspy tssCalling \
  -i sample1.bam -i sample2.bam -n "ctrl treat" \
  -r genome.fasta -o raw.TSS.tsv

# 2. Merge replicates + filter
tsspy mergeSamples merge      -i raw.TSS.tsv      -o merged.tsv -g "ctrl treat" -m "1 1 2 2"
tsspy mergeSamples filter     -i merged.tsv       -o filt.tsv   --method poisson --p-val 0.01 -r genome.fasta --normalization
# (or all-in-one)
tsspy mergeSamples process    -i raw.TSS.tsv      -o final.tsv  -g "ctrl treat" -m "1 1 2 2" --filter poisson -r genome.fasta

# 3. Cluster
tsspy clustering              -i filt.tsv         -o clusters --method peakclu
tsspy clustering              -i filt.tsv         -o clusters --method peakcluMax  # adamZhang fork extension

# 4. Cross-sample consensus
tsspy consensusCluster cluster -t filt.tsv \
  -i clusters.ctrl.tsv -i clusters.treat.tsv -n "ctrl treat" \
  -o consensus -d 50

# 5. Shape scores
tsspy shapeCluster batch \
  -c consensus.ctrl.tsv -c consensus.treat.tsv -n "ctrl treat" \
  -t filt.tsv -o shape -m PSS

# 6. Assign to genes
tsspy geneAssign assign \
  -c consensus.ctrl.tsv -c consensus.treat.tsv -n "ctrl treat" \
  -a annotation.gff -o assigned

# 7. Export to bedGraph / BigWig
tsspy bigwig -i filt.tsv --format bedGraph --batch --data-label processed
tsspy bigwig -i filt.tsv --format BigWig   --batch --data-label processed -r genome.fasta
```

## Two filter methods (matching TSSr)

```bash
# Poisson noise filter (needs raw counts; threshold derived from coverage)
tsspy mergeSamples filter -i merged.tsv -o filt.tsv \
  --method poisson --p-val 0.01 -r genome.fasta --normalization

# TPM threshold (needs normalized input)
tsspy mergeSamples filter -i normalized.tsv -o filt.tsv \
  --method TPM --tpm-low 0.1
```

## Differences from TSSr

| | TSSr (R) | TSSpy |
|---|---|---|
| BAM I/O | Rsamtools | pysam (read-by-read iterator + multiprocessing) |
| Reference genome | BSgenome (installed package per genome) | FASTA + `.fai` |
| BigWig export | rtracklayer | pyBigWig |
| CLI | R functions | `tsspy <command>` (typer-based) |
| Output values | matches TSSr exactly (in default mode) | matches TSSr bit-for-bit |

## Citation

If you use TSSpy, please cite TSSr alongside:
- TSSpy: https://github.com/JohnnyChen1113/TSSpy
- TSSr: Lu, Z., Berry, K., Hu, Z., Zhan, Y., Ahn, T., & Lin, Z. (2021). TSSr: an R package for comprehensive analyses of TSS sequencing data. *NAR Genomics and Bioinformatics*, 3(4).

## License

MIT
