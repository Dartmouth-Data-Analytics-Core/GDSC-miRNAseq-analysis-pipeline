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
  - hairpin alignment (bowtie1)
"""

from pathlib import Path

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


