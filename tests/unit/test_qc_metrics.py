"""
Unit tests for the log-parsing functions in scripts/qc_metrics-non-umi.py

The script file has a hyphen in its name so it cannot be imported normally.
We extract just the pure function definitions via AST and exec them.
BAM-dependent functions (run_samtools_count) are not tested here — they
require actual BAM files and are covered by the integration test.
"""

import ast
import textwrap
from pathlib import Path

import pytest

SCRIPTS = Path(__file__).parent.parent.parent / "scripts"
_SCRIPT = SCRIPTS / "qc_metrics-non-umi.py"


def _extract_functions(source: str, names: set) -> dict:
    """Exec only the specified function defs from source, return name→callable."""
    tree = ast.parse(source)
    ns: dict = {}
    for node in tree.body:
        if isinstance(node, ast.FunctionDef) and node.name in names:
            seg = ast.get_source_segment(source, node)
            exec(seg, ns)
    return ns


@pytest.fixture(scope="module")
def parse_fns():
    source = _SCRIPT.read_text()
    return _extract_functions(source, {"parse_bowtie2_log", "parse_bowtie1_log"})


def _write_bowtie2_log(path, total, aligned_once, multimapped):
    unaligned = total - aligned_once - multimapped
    path.write_text(textwrap.dedent(f"""\
        {total} reads; of these:
          {total} (100.00%) were unpaired; of these:
            {unaligned} ({unaligned/total*100:.2f}%) aligned 0 times
            {aligned_once} ({aligned_once/total*100:.2f}%) aligned exactly 1 time
            {multimapped} ({multimapped/total*100:.2f}%) aligned >1 times
        {(aligned_once+multimapped)/total*100:.2f}% overall alignment rate
    """))


def _write_bowtie1_log(path, total, mapped, unmapped):
    path.write_text(textwrap.dedent(f"""\
        # reads processed: {total}
        # reads with at least one alignment: {mapped}
        # reads that failed to align: {unmapped}
        Reported {mapped} alignments
    """))


class TestParseBowtie2Log:
    def test_correct_values(self, tmp_path, parse_fns):
        log = tmp_path / "sample.log"
        _write_bowtie2_log(log, total=10000, aligned_once=7800, multimapped=1000)
        total, mapped, unaligned, multimap = parse_fns["parse_bowtie2_log"](str(log))
        assert total == 10000
        assert mapped == 8800   # once + multi
        assert unaligned == 1200
        assert multimap == 1000

    def test_raises_on_empty_log(self, tmp_path, parse_fns):
        bad = tmp_path / "bad.log"
        bad.write_text("nothing useful\n")
        with pytest.raises(ValueError, match="Missing info"):
            parse_fns["parse_bowtie2_log"](str(bad))

    def test_all_unaligned(self, tmp_path, parse_fns):
        log = tmp_path / "zero.log"
        _write_bowtie2_log(log, total=5000, aligned_once=0, multimapped=0)
        total, mapped, unaligned, multimap = parse_fns["parse_bowtie2_log"](str(log))
        assert total == 5000
        assert mapped == 0
        assert unaligned == 5000


class TestParseBowtie1Log:
    def test_correct_values(self, tmp_path, parse_fns):
        log = tmp_path / "sample.log"
        _write_bowtie1_log(log, total=5000, mapped=3800, unmapped=1200)
        total, mapped, unmapped = parse_fns["parse_bowtie1_log"](str(log))
        assert total == 5000
        assert mapped == 3800
        assert unmapped == 1200

    def test_raises_on_missing_fields(self, tmp_path, parse_fns):
        bad = tmp_path / "bad.log"
        bad.write_text("no useful content\n")
        with pytest.raises(ValueError, match="Missing info"):
            parse_fns["parse_bowtie1_log"](str(bad))
