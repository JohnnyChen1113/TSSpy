"""End-to-end pipeline parity regression suite.

Each stage of the TSSpy pipeline is exercised on the 4-BAM yeast set and
diffed against the canonical TSSr 0.99.6 ground-truth artefacts in the
data dir. A failure here means we've regressed bit-for-bit (or, for
shape scores, FP-precision) parity with TSSr — which is the project's
stated north star (CLAUDE.md).
"""
from __future__ import annotations
import shutil
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from conftest import run_tsspy


# ---------- helpers ----------

def _read(path):
    return pd.read_csv(path, sep="\t")


def _rename_tssr(df):
    return df.rename(columns={"YPD.1": "s1", "YPD.2": "s2",
                              "Arrest.1": "s3", "Arrest.2": "s4"})


def _rename_tsspy(df):
    return df.rename(columns={"S01": "s1", "S02": "s2",
                              "S03": "s3", "S04": "s4"})


# ---------- Stage 1: tssCalling ----------

@pytest.fixture(scope="module")
def tssCalling_out(workdir, data_dir):
    out = workdir / "tsspy_stage1.TSS.tsv"
    run_tsspy(
        "tssCalling",
        "-i", str(data_dir / "S01.sorted.bam"),
        "-i", str(data_dir / "S02.sorted.bam"),
        "-i", str(data_dir / "S03.sorted.bam"),
        "-i", str(data_dir / "S04.sorted.bam"),
        "-n", "S01 S02 S03 S04",
        "-r", str(data_dir / "Scer.genome.fasta"),
        "-o", str(out),
    )
    return out


def test_tssCalling_bit_for_bit(tssCalling_out, data_dir):
    r = _rename_tssr(_read(data_dir / "tssr_repro_seqq10.TSS.tsv"))
    p = _rename_tsspy(_read(tssCalling_out))
    assert len(r) == len(p) == 163203
    m = r.merge(p, on=["chr", "pos", "strand"], how="outer", indicator=True)
    assert (m["_merge"] == "left_only").sum() == 0
    assert (m["_merge"] == "right_only").sum() == 0
    for col in ["s1", "s2", "s3", "s4"]:
        assert (r[col].values == p[col].values).all(), f"column {col} differs"


# ---------- Stage 2: mergeSamples merge ----------

@pytest.fixture(scope="module")
def merge_out(workdir, tssCalling_out, data_dir):
    out = workdir / "tsspy_stage2_merged.TSS.tsv"
    run_tsspy(
        "mergeSamples", "merge",
        "-i", str(tssCalling_out),
        "-o", str(out),
        "-g", "YPD Arrest",
        "-m", "1 1 2 2",
    )
    return out


def test_mergeSamples_bit_for_bit(merge_out, data_dir):
    r = _read(data_dir / "tssr_stage2_merged.TSS.tsv")
    p = _read(merge_out)
    assert len(r) == len(p) == 163203
    m = r.merge(p, on=["chr", "pos", "strand"], how="outer", suffixes=("_r", "_p"))
    for col in ["YPD", "Arrest"]:
        assert np.array_equal(m[f"{col}_r"].values, m[f"{col}_p"].values), f"column {col} differs"


# ---------- Stage 3a: normalize ----------

@pytest.fixture(scope="module")
def normalize_out(workdir, merge_out, data_dir):
    out = workdir / "tsspy_stage3a_normalized.TSS.tsv"
    run_tsspy(
        "mergeSamples", "normalize",
        "-i", str(merge_out),
        "-o", str(out),
    )
    return out


def test_normalize_bit_for_bit(normalize_out, data_dir):
    r = _read(data_dir / "tssr_stage3a_normalized.TSS.tsv")
    p = _read(normalize_out)
    assert len(r) == len(p)
    m = r.merge(p, on=["chr", "pos", "strand"], how="outer", suffixes=("_r", "_p"))
    for col in ["YPD", "Arrest"]:
        assert np.allclose(m[f"{col}_r"].values, m[f"{col}_p"].values, atol=1e-9), \
            f"column {col} differs beyond FP rounding"


# ---------- Stage 3b: poisson filter + normalize ----------

@pytest.fixture(scope="module")
def filter_out(workdir, merge_out, data_dir):
    out = workdir / "tsspy_stage3b_poisson_normalized.TSS.tsv"
    run_tsspy(
        "mergeSamples", "filter",
        "-i", str(merge_out),
        "-o", str(out),
        "--method", "poisson", "--p-val", "0.01",
        "--normalization",
        "-r", str(data_dir / "Scer.genome.fasta"),
    )
    return out


def test_filter_poisson_bit_for_bit(filter_out, data_dir):
    r = _read(data_dir / "tssr_stage3b_poisson_normalized.TSS.tsv")
    p = _read(filter_out)
    assert len(r) == len(p) == 116719
    m = r.merge(p, on=["chr", "pos", "strand"], how="outer", suffixes=("_r", "_p"))
    for col in ["YPD", "Arrest"]:
        assert np.array_equal(m[f"{col}_r"].values, m[f"{col}_p"].values), f"column {col} differs"


# ---------- Stage 4: clustering peakclu ----------

@pytest.fixture(scope="module")
def cluster_out(workdir, filter_out, data_dir):
    out_prefix = workdir / "tsspy_stage4_clusters"
    run_tsspy(
        "clustering",
        "-i", str(filter_out),
        "-o", str(out_prefix),
        "--method", "peakclu",
    )
    return out_prefix


@pytest.mark.parametrize("sample,expected_rows", [("YPD", 3307), ("Arrest", 3109)])
def test_clustering_peakclu_bit_for_bit(cluster_out, data_dir, sample, expected_rows):
    r = _read(data_dir / f"tssr_stage4_clusters_{sample}.tsv")
    p = _read(f"{cluster_out}.{sample}.tsv")
    assert len(r) == len(p) == expected_rows
    r_sorted = r.sort_values(["strand", "chr", "start"]).reset_index(drop=True)
    p_sorted = p.sort_values(["strand", "chr", "start"]).reset_index(drop=True)
    assert r_sorted.equals(p_sorted), f"{sample}: cluster DataFrames differ"


# ---------- Stage 5: consensusCluster ----------

@pytest.fixture(scope="module")
def consensus_out(workdir, filter_out, cluster_out, data_dir):
    out_prefix = workdir / "tsspy_stage5_consensus"
    run_tsspy(
        "consensusCluster", "cluster",
        "-t", str(filter_out),
        "-i", f"{cluster_out}.YPD.tsv",
        "-i", f"{cluster_out}.Arrest.tsv",
        "-n", "YPD Arrest",
        "-o", str(out_prefix),
        "-d", "50",
    )
    return out_prefix


@pytest.mark.parametrize("sample,expected_rows", [("YPD", 3306), ("Arrest", 3108)])
def test_consensus_bit_for_bit(consensus_out, data_dir, sample, expected_rows):
    r = _read(data_dir / f"tssr_stage5_consensus_{sample}.tsv")
    p = _read(f"{consensus_out}.{sample}.tsv")
    assert len(r) == len(p) == expected_rows
    r_sorted = r.sort_values(["strand", "chr", "start"]).reset_index(drop=True)
    p_sorted = p.sort_values(["strand", "chr", "start"]).reset_index(drop=True)
    assert r_sorted.equals(p_sorted)


# ---------- Stage 6: shapeCluster (PSS + SI) ----------

@pytest.fixture(scope="module")
def shape_out(workdir, filter_out, consensus_out, data_dir):
    for method in ("PSS", "SI"):
        out_prefix = workdir / f"tsspy_stage6_shape_{method}"
        run_tsspy(
            "shapeCluster", "batch",
            "-c", f"{consensus_out}.YPD.tsv",
            "-c", f"{consensus_out}.Arrest.tsv",
            "-n", "YPD Arrest",
            "-t", str(filter_out),
            "-o", str(out_prefix),
            "-m", method,
        )
    return workdir


@pytest.mark.parametrize("method,sample,expected_rows", [
    ("PSS", "YPD", 3306), ("PSS", "Arrest", 3108),
    ("SI",  "YPD", 3306), ("SI",  "Arrest", 3108),
])
def test_shape_fp_precision(shape_out, data_dir, method, sample, expected_rows):
    r = _read(data_dir / f"tssr_stage6_shape_{method}_{sample}.tsv")
    p = _read(shape_out / f"tsspy_stage6_shape_{method}.{sample}.tsv")
    assert len(r) == len(p) == expected_rows
    # Non-score columns must be exact
    for col in ["cluster", "chr", "start", "end", "strand", "dominant_tss",
                "tags", "tags.dominant_tss", "q_0.1", "q_0.9", "interquantile_width"]:
        if r[col].dtype == object:
            assert (r[col] == p[col]).all(), f"{method} {sample}: {col} differs"
        else:
            assert np.array_equal(r[col].values, p[col].values), f"{method} {sample}: {col} differs"
    # shape.score within FP precision (residual log/sum ordering noise R vs numpy)
    assert np.allclose(r["shape.score"].values, p["shape.score"].values, atol=1e-12)


# ---------- Stage 7: geneAssign ----------

@pytest.fixture(scope="module")
def assign_out(workdir, consensus_out, data_dir):
    out_prefix = workdir / "tsspy_stage7"
    run_tsspy(
        "geneAssign", "assign",
        "-c", f"{consensus_out}.YPD.tsv",
        "-c", f"{consensus_out}.Arrest.tsv",
        "-n", "YPD Arrest",
        "-a", str(data_dir / "saccharomyces_cerevisiae_R64-2-1.gff"),
        "-o", str(out_prefix),
        "--filter-cluster", "--filter-cluster-threshold", "0.02",
    )
    return out_prefix


@pytest.mark.parametrize("kind,sample,expected_rows", [
    ("assigned",   "YPD",    860),
    ("assigned",   "Arrest", 766),
    ("unassigned", "YPD",    2446),
    ("unassigned", "Arrest", 2342),
])
def test_geneAssign_bit_for_bit(assign_out, data_dir, kind, sample, expected_rows):
    r = _read(data_dir / f"tssr_stage7_{kind}_{sample}.tsv")
    p = _read(f"{assign_out}.{sample}.{kind}.tsv")
    assert len(r) == len(p) == expected_rows
    # cluster-set match
    assert set(r["cluster"]) == set(p["cluster"])
    # All common columns identical (string-compared so NA handled uniformly)
    r_s = r.sort_values("cluster").reset_index(drop=True)
    p_s = p.sort_values("cluster").reset_index(drop=True)
    common = [c for c in r_s.columns if c in p_s.columns]
    for c in common:
        assert (r_s[c].fillna("__NA__").astype(str)
                == p_s[c].fillna("__NA__").astype(str)).all(), \
            f"{kind} {sample}: column {c} differs"


# ---------- Stage 8: bigwig export ----------

@pytest.fixture(scope="module")
def bedgraph_out(workdir, filter_out, data_dir):
    out_dir = workdir / "bw"
    out_dir.mkdir(exist_ok=True)
    # The CLI's --batch mode writes to CWD; copy filter input then cd in subprocess
    # Using shutil to drive directly via the module API would be cleaner; subprocess
    # with cwd= is the simplest portable equivalent.
    run_tsspy(
        "bigwig",
        "-i", str(filter_out),
        "--format", "bedGraph",
        "--batch", "--data-label", "processed",
        cwd=str(out_dir),
    )
    return out_dir


@pytest.mark.parametrize("sample,strand", [
    ("YPD", "plus"), ("YPD", "minus"),
    ("Arrest", "plus"), ("Arrest", "minus"),
])
def test_bedgraph_byte_for_byte(bedgraph_out, data_dir, sample, strand):
    """bedGraph TSV is a plain text format → expect byte-identical files."""
    r_path = data_dir / "tssr_stage8_bw" / f"{sample}.TSS.processed.{strand}.bedGraph"
    p_path = bedgraph_out / f"{sample}.TSS.processed.{strand}.bedGraph"
    r_bytes = r_path.read_bytes()
    p_bytes = p_path.read_bytes()
    assert r_bytes == p_bytes, f"{sample} {strand}: bedGraph differs by {len(r_bytes)-len(p_bytes):+d} bytes"
