"""
Unit tests for scripts/add_biotype.py

Tests that:
- gene_biotype and gene_name columns are appended correctly
- genes missing from GTF get empty string biotype/name (not NaN crash)
- column ordering is preserved (Geneid, gene_name, gene_biotype, then rest)
- output file is written correctly
"""

import subprocess
import sys
import textwrap
from pathlib import Path

import pandas as pd
import pytest

SCRIPT = Path(__file__).parent.parent.parent / "scripts" / "add_biotype.py"


def run_script(fc_path, gtf_path, out_path):
    result = subprocess.run(
        [sys.executable, str(SCRIPT), str(fc_path), str(gtf_path), str(out_path)],
        capture_output=True,
        text=True,
    )
    return result


def test_biotype_added(tmp_path, sample_featurecounts_tsv, sample_gtf):
    out = tmp_path / "out.tsv"
    result = run_script(sample_featurecounts_tsv, sample_gtf, out)
    assert result.returncode == 0, result.stderr
    df = pd.read_csv(out, sep="\t")
    assert "gene_biotype" in df.columns
    assert "gene_name" in df.columns


def test_column_order(tmp_path, sample_featurecounts_tsv, sample_gtf):
    out = tmp_path / "out.tsv"
    run_script(sample_featurecounts_tsv, sample_gtf, out)
    df = pd.read_csv(out, sep="\t")
    cols = df.columns.tolist()
    assert cols[0] == "Geneid"
    assert cols[1] == "gene_name"
    assert cols[2] == "gene_biotype"


def test_known_biotypes(tmp_path, sample_featurecounts_tsv, sample_gtf):
    out = tmp_path / "out.tsv"
    run_script(sample_featurecounts_tsv, sample_gtf, out)
    df = pd.read_csv(out, sep="\t")
    row_mirna = df[df["Geneid"] == "ENSDARG00000001"]
    assert row_mirna["gene_biotype"].values[0] == "miRNA"
    assert row_mirna["gene_name"].values[0] == "dre-mir-21"


def test_missing_gene_gets_empty_string(tmp_path, sample_featurecounts_tsv, sample_gtf):
    """Gene in featureCounts but absent from GTF should get empty string, not NaN crash."""
    extra_fc = textwrap.dedent("""\
        # Program:featureCounts
        Geneid\tChr\tStart\tEnd\tStrand\tLength\tsample_1.bam
        ENSDARG99999999\tchr1\t5000\t6000\t+\t1000\t50
    """)
    fc_path = tmp_path / "extra_fc.tsv"
    fc_path.write_text(extra_fc)
    out = tmp_path / "out.tsv"
    result = run_script(fc_path, sample_gtf, out)
    assert result.returncode == 0
    df = pd.read_csv(out, sep="\t")
    row = df[df["Geneid"] == "ENSDARG99999999"]
    assert row["gene_biotype"].values[0] == "" or pd.isna(row["gene_biotype"].values[0])


def test_output_row_count(tmp_path, sample_featurecounts_tsv, sample_gtf):
    out = tmp_path / "out.tsv"
    run_script(sample_featurecounts_tsv, sample_gtf, out)
    df = pd.read_csv(out, sep="\t")
    assert len(df) == 2
