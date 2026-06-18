"""
Integration tests: Pipeline output validation

Run AFTER snakemake completes in the integration CI workflow.
Checks that expected output files exist, are non-empty, and have
the correct format/columns.

These tests assume the pipeline was run from the repository root
with tests/fixtures/test_config.yaml targeting:
  - trimming
  - FastQC
  - seqcluster
  - mirtop (isomiR path)
  - pivot + collation (R scripts)
"""

from pathlib import Path

import pandas as pd
import pytest

REPO = Path(__file__).parent.parent.parent
SAMPLES = ["sample_1", "sample_2"]


# ── Helpers ───────────────────────────────────────────────────────────────────

def assert_file_exists_nonempty(path):
    p = REPO / path
    assert p.exists(), f"Expected output missing: {path}"
    assert p.stat().st_size > 0, f"Output file is empty: {path}"


# ── Trimming ──────────────────────────────────────────────────────────────────

class TestTrimmingOutputs:
    @pytest.mark.parametrize("sample", SAMPLES)
    def test_trimmed_fastq_exists(self, sample):
        assert_file_exists_nonempty(f"trimming/{sample}.R1.trim.fastq.gz")

    @pytest.mark.parametrize("sample", SAMPLES)
    def test_cutadapt_log_exists(self, sample):
        assert_file_exists_nonempty(f"trimming/logs/{sample}.cutadapt.log")


# ── FastQC ────────────────────────────────────────────────────────────────────

class TestFastQCOutputs:
    def test_multiqc_config_created(self):
        assert_file_exists_nonempty("fastQC/fastqc_multiqc_config.yaml")

    @pytest.mark.parametrize("sample", SAMPLES)
    def test_fastqc_html_exists(self, sample):
        assert_file_exists_nonempty(f"fastQC/{sample}.R1.trim_fastqc.html")

    @pytest.mark.parametrize("sample", SAMPLES)
    def test_fastqc_zip_exists(self, sample):
        assert_file_exists_nonempty(f"fastQC/{sample}.R1.trim_fastqc.zip")


# ── SeqCluster ────────────────────────────────────────────────────────────────

class TestSeqclusterOutputs:
    @pytest.mark.parametrize("sample", SAMPLES)
    def test_collapsed_fastq_exists(self, sample):
        assert_file_exists_nonempty(f"collapsed/{sample}.seqcluster.fastq.gz")


# ── Hairpin alignment ─────────────────────────────────────────────────────────

class TestHairpinAlignmentOutputs:
    @pytest.mark.parametrize("sample", SAMPLES)
    def test_hairpin_bam_exists(self, sample):
        assert_file_exists_nonempty(f"collapsed/{sample}.seqcluster.hairpin.aln.srt.bam")

    @pytest.mark.parametrize("sample", SAMPLES)
    def test_hairpin_idxstats_exists(self, sample):
        assert_file_exists_nonempty(f"collapsed/{sample}.seqcluster.hairpin.aln.srt.bam.idxstats")

    @pytest.mark.parametrize("sample", SAMPLES)
    def test_hairpin_flagstat_exists(self, sample):
        assert_file_exists_nonempty(f"collapsed/{sample}.seqcluster.hairpin.aln.srt.bam.flagstat")


# ── miRtop ────────────────────────────────────────────────────────────────────

class TestMirtopOutputs:
    @pytest.mark.parametrize("sample", SAMPLES)
    def test_mirtop_gff_exists(self, sample):
        assert_file_exists_nonempty(f"mirtop/{sample}.hairpin.gff")

    @pytest.mark.parametrize("sample", SAMPLES)
    def test_mirtop_tsv_exists(self, sample):
        assert_file_exists_nonempty(f"mirtop/{sample}.hairpin.tsv")

    def test_mirtop_stats_log_exists(self):
        assert_file_exists_nonempty("mirtop/mirtop_stats.log")


# ── Count table (isomiR path) ─────────────────────────────────────────────────

class TestIsomirCountTable:
    def test_master_table_exists(self):
        assert_file_exists_nonempty("miRNA_Quant/raw_merged_canonical_and_all_isomirs.csv")

    def test_master_table_has_rows(self):
        df = pd.read_csv(REPO / "miRNA_Quant/raw_merged_canonical_and_all_isomirs.csv")
        assert len(df) > 0, "isomiR count table is empty"

    def test_master_table_has_sample_columns(self):
        df = pd.read_csv(REPO / "miRNA_Quant/raw_merged_canonical_and_all_isomirs.csv")
        for sample in SAMPLES:
            assert sample in df.columns, f"Sample column {sample} missing from count table"

    def test_master_table_has_mirna_column(self):
        df = pd.read_csv(REPO / "miRNA_Quant/raw_merged_canonical_and_all_isomirs.csv")
        assert "miRNA" in df.columns, "miRNA column missing from count table"

    def test_canonical_counts_file_exists(self):
        assert_file_exists_nonempty("miRNA_Quant/raw_canonical_counts.tsv")
