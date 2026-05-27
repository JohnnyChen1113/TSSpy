---
name: tsspy
description: Use when working with TSSpy (the Python port of TSSr) — pipeline help, CLI usage, parity guarantees, the --strict-tssr flag, and three known TSSr bugs documented for upstream. Triggers on mentions of TSSpy, TSSr, CAGE-seq analysis, TSS calling, peakclu / consensusCluster / shapeCluster / annotateCluster / exportTSStoBedgraph, or when working in a directory containing BAM + TSSr-style outputs.
---

# TSSpy — Python port of TSSr at bit-for-bit parity

TSSpy is a Python CLI that reproduces TSSr 0.99.6's full TSS-analysis pipeline. The default mode (`--strict-tssr=true`) produces output byte-identical to TSSr on the same input; the opt-in `--no-strict-tssr` mode fixes two known TSSr bugs.

## Pipeline (7 stages, all at TSSr parity)

```
tssCalling  →  mergeSamples merge  →  normalize  →  filter  →
clustering  →  consensusCluster  →  shapeCluster  →  geneAssign  →  bigwig
```

Plus diagnostic plots:
```
correlation  |  plot pca / iqw / shape / tss
```

## CLI quick reference

```bash
# 1. Call TSSs from BAMs (default: --strict-tssr=true, sequencing-quality=10)
tsspy tssCalling \
  -i s1.bam -i s2.bam -n "ctrl treat" \
  -r genome.fasta -o raw.TSS.tsv

# 2. Merge replicates
tsspy mergeSamples merge -i raw.TSS.tsv -o merged.tsv \
  -g "ctrl treat" -m "1 1 2 2"

# 3. Filter + normalize (poisson noise filter, TSSr default)
tsspy mergeSamples filter -i merged.tsv -o filt.tsv \
  --method poisson --p-val 0.01 -r genome.fasta --normalization
# Alternative: TPM-based filter (requires normalized input)
tsspy mergeSamples filter -i normalized.tsv -o filt.tsv \
  --method TPM --tpm-low 0.1

# 4. Cluster
tsspy clustering -i filt.tsv -o clusters --method peakclu
# Or use Adam Zhang's peakcluMax variant
tsspy clustering -i filt.tsv -o clusters --method peakcluMax

# 5. Cross-sample consensus
tsspy consensusCluster cluster \
  -t filt.tsv \
  -i clusters.ctrl.tsv -i clusters.treat.tsv -n "ctrl treat" \
  -o consensus -d 50

# 6. Promoter shape scores (PSS or SI)
tsspy shapeCluster batch \
  -c consensus.ctrl.tsv -c consensus.treat.tsv -n "ctrl treat" \
  -t filt.tsv -o shape -m PSS

# 7. Assign to genes (GFF3)
tsspy geneAssign assign \
  -c consensus.ctrl.tsv -c consensus.treat.tsv -n "ctrl treat" \
  -a annotation.gff -o assigned

# 8. Export to bedGraph / BigWig
tsspy bigwig -i filt.tsv --format bedGraph --batch --data-label processed
tsspy bigwig -i filt.tsv --format BigWig   --batch --data-label processed -r genome.fasta

# Diagnostic plots (plotnine-based, ggplot2-style)
tsspy plot pca   -i raw.TSS.tsv -o pca.pdf --tss-threshold 10 \
                 --merge-labels "ctrl treat" --merge-index "1 1 2 2"
tsspy plot iqw   -i consensus.ctrl.tsv -i consensus.treat.tsv \
                 -n "ctrl treat" -o iqw.pdf
tsspy plot shape -i shape.ctrl.tsv -i shape.treat.tsv \
                 -n "ctrl treat" -o shape.pdf
tsspy plot tss   -t filt.tsv -c consensus.ctrl.tsv -c consensus.treat.tsv \
                 -n "ctrl treat" -a annotation.gff \
                 --genes "YAL003W YBR090C" -o tss_browser.pdf
tsspy correlation -i raw.TSS.tsv --plot --plot-file corr.png --source raw
```

## Default parameter values (all match TSSr 0.99.6)

| Param | Default | Where |
|---|---|---|
| `--sequencing-quality` | 10 | tssCalling |
| `--mapping-quality` | 20 | tssCalling |
| `--strict-tssr` | true | tssCalling |
| `--p-val` (poisson) | 0.01 | filter |
| `--tpm-low` | 0.1 | filter |
| `--peak-distance` | 100 | clustering |
| `--extension-distance` | 30 | clustering |
| `--local-threshold` | 0.02 | clustering |
| `--cluster-threshold` | 1.0 | clustering |
| `-d` / `--dis` | 50 | consensusCluster |
| `--upstream` | 1000 | geneAssign |
| `--upstream-overlap` | 500 | geneAssign |
| `--downstream` | 0 | geneAssign |
| `--filter-cluster-threshold` | 0.02 | geneAssign |

## --strict-tssr=true vs false — what each mode does

**Default (`true`)**: reproduces TSSr 0.99.6 byte-for-byte. Preserves three known TSSr behaviours:

1. **mapped_length sums all CIGAR ints** (incl. I and S). For insertion-containing minus-strand reads this shifts the TSS 1bp past the actual alignment end.
2. **Minus-strand G-mismatch iteration checks the wrong end of the read** (cDNA 3' end via `seq[1..i]` instead of cDNA 5' end via `seq[last-i+1..last]`).
3. Plus-strand bug-free.

**Corrected (`--no-strict-tssr`)**: applies the two fixes. Use when you want biologically faithful TSS positions but accept divergence from TSSr ground truth.

**Validation evidence (in TSSr_fix_proposal.md)**:
- S. cerevisiae S288C (4 BAM): read-weighted YR Inr-motif rate 84.06% → 85.14% (+1.08 pp), minus-strand only 82.69% → 84.77% (+2.08 pp).
- S. uvarum (2 BAM, independent species): 52.05% → 53.79% (+1.74 pp), minus 49.91% → 53.64% (+3.72 pp). Replicates the effect.

## Three TSSr bugs documented for upstream

Full write-up: `TSSr_fix_proposal.md` (repo root).

| # | TSSr location | Symptom |
|---|---|---|
| 1 | `R/ImportFunctions.R:62-69` (`getTSS`) | mapped_length overcounts I → minus-strand TSS off by 1bp on ~6,857 insertion reads (S288C) |
| 2 | `R/ImportFunctions.R:155` (`.removeNewG`) | Minus-strand iteration filter on wrong end of read → ~104K reads over/under-corrected |
| 3 | `R/AnnotationMethods.R:128` (`annotateCluster`) | `rbind(m[,seq(12)], ...)` data.table j-expression trap → `@filteredClusters` output replaced by 12-column garbage matrix |

## Output file format conventions

All TSS-table files: tab-separated, columns `chr pos strand <sample1> [<sample2> ...]`. Positions are 1-based.

All cluster files (peakclu / consensusCluster output): `cluster chr start end strand dominant_tss tags tags.dominant_tss q_0.1 q_0.9 interquantile_width`.

Shape output appends `shape.score` to the cluster file.

geneAssign output appends `gene` (and optionally `inCoding`) to the cluster file. Split into `<prefix>.<sample>.assigned.tsv` and `.unassigned.tsv`.

bedGraph: `chr  pos-1  pos  value` (no track header, minus-strand values negated). One file per sample × strand named `<sample>.TSS.<data>.{plus,minus}.bedGraph`.

## Verifying parity

The canonical ground-truth file for `tssCalling` is `tssr_repro_seqq10.TSS.tsv`. Regenerate by running `Rscript repro_tssr.R` from repo root (~3 min). Per-stage ground-truth files live alongside (`tssr_stage{2,3,4,5,6,7,8}_*.tsv`); regenerate by running `Rscript repro_tssr_downstream.R` etc.

**Do NOT use `ALL.samples.TSS.raw.txt_github` as ground truth** — that file is not reproducible from TSSr 0.99.6 with default params (~5,000 positions off; origin unknown).

Run `pytest tests/` from the repo root to verify the full 24-test parity suite (~94s).

## Common gotchas

- **Sample column names must match between TSS table and cluster files.** consensusCluster, shapeCluster, plot iqw/shape all expect the per-sample cluster TSV to share its sample name with a column in the TSS table.
- **GFF gene IDs may be URL-encoded** (e.g. `tP%28UGG%29A`). TSSpy decodes via `urllib.parse.unquote`. Hand-curated annotations with literal `(` / `)` work too.
- **TSSpy expects gene-like feature types in GFF**: `gene`, `tRNA_gene`, `snoRNA_gene`, `rRNA_gene`, `ncRNA_gene`, `pseudogene`, `snRNA_gene`, `telomerase_RNA_gene`, `transposable_element_gene`.
- **shape.score has FP-precision discrepancy with TSSr** of up to ~1e-12 due to log/sum ordering between R and numpy. All other columns are bit-for-bit.
- **TSSr's `@filteredClusters` is broken upstream**; TSSpy emits the correct output (~2,150 / 2,096 filtered clusters for YPD / Arrest). Direct diff against TSSr's filtered output will not work.

## When asked to add or modify pipeline behavior

1. **Read TSSr's R source first** (in `TSSr/R/`) — TSSpy mirrors TSSr exactly by design.
2. **Check `TSSr_fix_proposal.md`** for known bugs we deliberately preserve in `--strict-tssr` mode.
3. **After any change, run `pytest tests/`** to confirm no parity regression.
4. **For new analyses**, mirror TSSr's algorithm exactly in `--strict-tssr=true`; if introducing a bug-fix mode, add a flag and document the deviation.

## Architecture notes

- Each pipeline module exposes `app = typer.Typer()` and is registered in `TSSpy/main.py`.
- All quantile computations use integer-scaled arithmetic (tags × 1e6 → int64) to avoid FP ties at the 10% boundary.
- pandas `sort_values(kind='stable')` is required to match R's `data.table::setorder` semantics when sort keys have ties (used in `gene_assign._promoter_regions` and elsewhere).
- Multiprocessing in tssCalling uses one process per BAM, matching TSSr's `mclapply` behaviour.

## Repository

- Code: https://github.com/JohnnyChen1113/TSSpy
- TSSr upstream: https://github.com/Linlab-slu/TSSr
- adamZhang_TSSr fork (peakcluMax + faster BAM reader): `adamZhang_TSSr/` colocated in repo
