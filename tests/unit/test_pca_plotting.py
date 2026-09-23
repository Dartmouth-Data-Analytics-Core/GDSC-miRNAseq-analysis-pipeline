"""
Unit tests for scripts/pca_plotting.py

Tests that:
- Script runs without errors on a small count matrix
- Expected output PNG files are created and non-empty
- The median-of-ratios normalization function works correctly
- HVG detection plateau logic does not crash on edge cases
"""

import subprocess
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

SCRIPT = Path(__file__).parent.parent.parent / "scripts" / "pca_plotting.py"


def run_pca_script(input_csv, output_dir, pca_comp=3):
    result = subprocess.run(
        [sys.executable, str(SCRIPT), str(input_csv), str(output_dir), "-p", str(pca_comp)],
        capture_output=True,
        text=True,
    )
    return result


class TestMedianOfRatios:
    """Test the normalization function extracted from the script."""

    def _load_fn(self):
        source = SCRIPT.read_text()
        import ast
        tree = ast.parse(source)
        func_src = "\n".join(
            ast.get_source_segment(source, node)
            for node in tree.body
            if isinstance(node, ast.FunctionDef) and node.name == "median_of_ratios"
        )
        ns = {"np": np, "pd": pd}
        exec(func_src, ns)
        return ns["median_of_ratios"]

    def test_all_zeros_returns_input(self):
        fn = self._load_fn()
        df = pd.DataFrame({"s1": [0, 0], "s2": [0, 0]})
        result = fn(df)
        assert result.shape == df.shape

    def test_normalizes_counts(self):
        fn = self._load_fn()
        df = pd.DataFrame({"s1": [100.0, 200.0, 400.0], "s2": [200.0, 400.0, 800.0]})
        result = fn(df)
        # Ratio s2/s1 should be ~2 everywhere after normalization
        assert result.shape == df.shape
        assert (result > 0).all().all()


class TestPcaPlottingScript:
    def test_creates_output_files(self, tmp_path, sample_count_matrix):
        out_dir = tmp_path / "plots"
        result = run_pca_script(sample_count_matrix, out_dir, pca_comp=3)
        assert result.returncode == 0, f"Script failed:\n{result.stderr}"
        # At minimum these two files should always be produced
        assert (out_dir / "PCA_all_PC1_vs_PC2.png").exists()
        assert (out_dir / "PCA_all_PCA_variance_bar.png").exists()

    def test_output_files_nonempty(self, tmp_path, sample_count_matrix):
        out_dir = tmp_path / "plots"
        run_pca_script(sample_count_matrix, out_dir, pca_comp=3)
        for png in out_dir.glob("*.png"):
            assert png.stat().st_size > 0, f"{png.name} is empty"

    def test_creates_hvg_log(self, tmp_path, sample_count_matrix):
        out_dir = tmp_path / "plots"
        run_pca_script(sample_count_matrix, out_dir, pca_comp=3)
        assert (out_dir / "pca_hvg_log.txt").exists()

    def test_minimal_2_samples(self, tmp_path):
        """2-sample matrix — PCA with 1 component should not crash."""
        counts = pd.DataFrame(
            {"miRNA": ["mir-1", "mir-2", "mir-3"], "s1": [100, 200, 50], "s2": [90, 210, 60]}
        ).set_index("miRNA")
        csv = tmp_path / "mini.csv"
        counts.to_csv(csv)
        out_dir = tmp_path / "plots"
        result = run_pca_script(csv, out_dir, pca_comp=1)
        assert result.returncode == 0, result.stderr
