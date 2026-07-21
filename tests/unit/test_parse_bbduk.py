"""
Unit tests for scripts/parse_bbduk_spikein_stats.py

Tests the parse_bbduk_stats() function and the full main() pipeline.
"""

import subprocess
import sys
import textwrap
from pathlib import Path

import pandas as pd
import pytest

SCRIPT = Path(__file__).parent.parent.parent / "scripts" / "parse_bbduk_spikein_stats.py"


def _load_parse_fn():
    source = SCRIPT.read_text()
    import ast
    tree = ast.parse(source)
    # Collect everything needed: EXPECTED_SPIKEINS constant + parse_bbduk_stats fn
    needed_names = {"EXPECTED_SPIKEINS", "parse_bbduk_stats"}
    pieces = []
    for node in tree.body:
        if isinstance(node, ast.Assign):
            for target in node.targets:
                if isinstance(target, ast.Name) and target.id in needed_names:
                    pieces.append(ast.get_source_segment(source, node))
        elif isinstance(node, ast.FunctionDef) and node.name in needed_names:
            pieces.append(ast.get_source_segment(source, node))
    ns = {"sys": sys, "re": __import__("re")}
    exec("\n".join(pieces), ns)
    return ns["parse_bbduk_stats"], ns["EXPECTED_SPIKEINS"]


def make_stats_file(tmp_path, name, counts):
    """Write a BBDuk-format stats file."""
    lines = ["#File\tsample.fastq.gz", "#Total\t100000", "#Matched\t10000",
             "#Name\tReads\tReadsPct"]
    for spike_id, n in counts.items():
        lines.append(f"#{spike_id}\t{n}\t0.01%")
    f = tmp_path / name
    f.write_text("\n".join(lines) + "\n")
    return f


class TestParseBbdukStats:
    def test_all_spikeins_present(self, tmp_path):
        parse_fn, EXPECTED = _load_parse_fn()
        counts = {f"miND-{i:02d}": i * 100 for i in range(1, 8)}
        f = make_stats_file(tmp_path, "s1.stats", counts)
        result = parse_fn(str(f))
        assert set(result.keys()) == set(EXPECTED)
        for i in range(1, 8):
            assert result[f"miND-{i:02d}"] == i * 100

    def test_missing_spikein_filled_zero(self, tmp_path):
        parse_fn, EXPECTED = _load_parse_fn()
        counts = {"miND-01": 500, "miND-03": 300}  # miND-02 missing
        f = make_stats_file(tmp_path, "s1.stats", counts)
        result = parse_fn(str(f))
        assert result["miND-02"] == 0

    def test_empty_file_gives_zeros(self, tmp_path):
        parse_fn, EXPECTED = _load_parse_fn()
        f = tmp_path / "empty.stats"
        f.write_text("#File\t-\n#Total\t0\n#Matched\t0\n")
        result = parse_fn(str(f))
        for spike_id in EXPECTED:
            assert result[spike_id] == 0


class TestParseBbdukMain:
    def _run(self, tmp_path, stats_files, samples, output_name="out.tsv"):
        out = tmp_path / output_name
        cmd = [
            sys.executable, str(SCRIPT),
            "--stats", *[str(f) for f in stats_files],
            "--samples", *samples,
            "--output", str(out),
        ]
        result = subprocess.run(cmd, capture_output=True, text=True)
        return result, out

    def test_two_samples_output_shape(self, tmp_path):
        s1 = make_stats_file(tmp_path, "s1.stats", {"miND-01": 100, "miND-02": 200})
        s2 = make_stats_file(tmp_path, "s2.stats", {"miND-01": 150, "miND-02": 250})
        result, out = self._run(tmp_path, [s1, s2], ["sample_1", "sample_2"])
        assert result.returncode == 0, result.stderr
        df = pd.read_csv(out, sep="\t", index_col=0)
        assert df.shape == (7, 2)  # 7 spike-ins × 2 samples
        assert "sample_1" in df.columns
        assert df.loc["miND-01", "sample_1"] == 100

    def test_mismatched_counts_raises(self, tmp_path):
        s1 = make_stats_file(tmp_path, "s1.stats", {})
        result, _ = self._run(tmp_path, [s1], ["sample_1", "sample_2"])
        assert result.returncode != 0
