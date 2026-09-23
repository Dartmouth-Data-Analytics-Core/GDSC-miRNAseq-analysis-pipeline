"""
Unit tests for scripts/readcnt_to_rpkmtpm.py

Tests RPKM and TPM normalization functions, and that the script writes
the correct output files from a featureCounts-style TSV.
"""

import subprocess
import sys
import textwrap
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

SCRIPT = Path(__file__).parent.parent.parent / "scripts" / "readcnt_to_rpkmtpm.py"


def _load_fns():
    source = SCRIPT.read_text()
    import ast
    tree = ast.parse(source)
    ns = {"np": np, "pd": pd, "sys": sys}
    for node in tree.body:
        if isinstance(node, ast.FunctionDef) and node.name in {"to_rpkm_pre_sliced", "to_tpm_pre_sliced"}:
            exec(ast.get_source_segment(source, node), ns)
    return ns["to_rpkm_pre_sliced"], ns["to_tpm_pre_sliced"]


def make_fc_biotype_tsv(tmp_path, n_genes=5, n_samples=3):
    """Write a featureCounts output TSV with biotype columns as expected by the script."""
    rng = np.random.default_rng(0)
    counts = rng.integers(50, 2000, size=(n_genes, n_samples))
    lengths = rng.integers(500, 5000, size=n_genes)

    header = ["Geneid", "gene_name", "gene_biotype", "Chr", "Start", "End", "Strand", "Length"]
    sample_cols = [f"genome_alignment/s{i}.genome.srt.filt.bam" for i in range(n_samples)]
    header += sample_cols

    rows = []
    for i in range(n_genes):
        row = [f"ENSG{i:07d}", f"gene_{i}", "miRNA", "chr1",
               str(i * 1000 + 1), str(i * 1000 + lengths[i]), "+", str(lengths[i])]
        row += [str(c) for c in counts[i]]
        rows.append("\t".join(row))

    content = "\t".join(header) + "\n" + "\n".join(rows) + "\n"
    f = tmp_path / "featurecounts.readcounts.biotype.tsv"
    f.write_text(content)
    return f, lengths, counts


class TestNormalizationFunctions:
    def test_rpkm_shape(self):
        fn_rpkm, _ = _load_fns()
        arr = np.array([[1000.0, 100.0, 200.0],
                        [2000.0, 150.0, 300.0]])
        result = fn_rpkm(arr)
        assert result.shape == (2, 2)  # 2 genes × 2 samples

    def test_tpm_sums_to_million(self):
        _, fn_tpm = _load_fns()
        arr = np.array([[1000.0, 100.0, 200.0],
                        [2000.0, 150.0, 300.0],
                        [500.0, 50.0, 100.0]])
        result = fn_tpm(arr)
        col_sums = result.sum(axis=0)
        np.testing.assert_allclose(col_sums, 1e6, rtol=1e-4)

    def test_rpkm_longer_gene_lower_value(self):
        fn_rpkm, _ = _load_fns()
        # 2 genes, same count but different lengths
        arr = np.array([[1000.0, 100.0],   # short gene
                        [10000.0, 100.0]])  # 10× longer gene, same count
        result = fn_rpkm(arr)
        assert result[0, 0] > result[1, 0]


class TestRpkmScript:
    def _run(self, tmp_path, tsv_path, layout="single"):
        result = subprocess.run(
            [sys.executable, str(SCRIPT), str(tsv_path), layout],
            capture_output=True,
            text=True,
        )
        return result

    def test_creates_rpkm_and_tpm_files(self, tmp_path):
        tsv, _, _ = make_fc_biotype_tsv(tmp_path)
        result = self._run(tmp_path, tsv)
        assert result.returncode == 0, result.stderr
        stem = str(tsv)[:-4]
        assert Path(f"{stem}_rpkm.tsv").exists()
        assert Path(f"{stem}_tpm.tsv").exists()

    def test_paired_creates_fpkm_file(self, tmp_path):
        tsv, _, _ = make_fc_biotype_tsv(tmp_path)
        result = self._run(tmp_path, tsv, layout="paired")
        assert result.returncode == 0, result.stderr
        stem = str(tsv)[:-4]
        assert Path(f"{stem}_fpkm.tsv").exists()

    def test_output_annotation_columns_preserved(self, tmp_path):
        tsv, _, _ = make_fc_biotype_tsv(tmp_path)
        self._run(tmp_path, tsv)
        stem = str(tsv)[:-4]
        df = pd.read_csv(f"{stem}_rpkm.tsv", sep="\t")
        assert "Geneid" in df.columns
        assert "gene_biotype" in df.columns
