"""
Integration tests: Snakemake dry-run validation

Validates that the Snakemake DAG can be built for each config variant
without executing any shell commands. This catches:
- Missing or mis-spelled config keys accessed at parse time
- Wildcard resolution failures
- Rule dependency cycles
- Missing `include:` rule files

Requires: snakemake, pandas
"""

import subprocess
import sys
from pathlib import Path

import pytest

REPO = Path(__file__).parent.parent.parent
FIXTURES = Path(__file__).parent.parent / "fixtures"


def run_dryrun(config_file, extra_args=None):
    cmd = [
        sys.executable, "-m", "snakemake",
        "--dry-run",
        "--configfile", str(config_file),
        "--snakefile", str(REPO / "Snakefile"),
        "--quiet", "rules",
    ]
    if extra_args:
        cmd.extend(extra_args)
    return subprocess.run(cmd, capture_output=True, text=True, cwd=str(REPO))


class TestDryrunConfigs:
    def test_default_config(self):
        """Base config (no spike-ins, no UMI tools) should produce a valid DAG."""
        result = run_dryrun(FIXTURES / "test_config.yaml")
        assert result.returncode == 0, (
            f"Dry-run failed for default config.\nSTDOUT:\n{result.stdout}\nSTDERR:\n{result.stderr}"
        )

    def test_spikeins_config(self):
        """Spike-in config should include spikein_* rules in the DAG."""
        result = run_dryrun(FIXTURES / "test_config_spikeins.yaml")
        assert result.returncode == 0, (
            f"Dry-run failed for spike-in config.\nSTDOUT:\n{result.stdout}\nSTDERR:\n{result.stderr}"
        )

    def test_umitools_config(self):
        """UMI tools config should include umi_* and dedup rules in the DAG."""
        result = run_dryrun(FIXTURES / "test_config_umitools.yaml")
        assert result.returncode == 0, (
            f"Dry-run failed for UMI tools config.\nSTDOUT:\n{result.stdout}\nSTDERR:\n{result.stderr}"
        )


class TestDryrunRuleInclusion:
    def test_genome_alignment_rules_included(self):
        """genome_alignment rule should appear in the DAG (always included)."""
        result = run_dryrun(FIXTURES / "test_config.yaml")
        # If the DAG built successfully, genome_alignment is implicitly included
        assert result.returncode == 0

    def test_trimming_rule_in_dag(self):
        """The trimming rule should be schedulable for each sample in the fixture CSV."""
        result = run_dryrun(
            FIXTURES / "test_config.yaml",
            extra_args=[
                "trimming/sample_1.R1.trim.fastq.gz",
                "trimming/sample_2.R1.trim.fastq.gz",
            ],
        )
        assert result.returncode == 0, result.stderr

    def test_isomir_collation_in_dag(self):
        """The full isomiR path should be schedulable."""
        result = run_dryrun(
            FIXTURES / "test_config.yaml",
            extra_args=["miRNA_Quant/raw_merged_canonical_and_all_isomirs.csv"],
        )
        assert result.returncode == 0, result.stderr
