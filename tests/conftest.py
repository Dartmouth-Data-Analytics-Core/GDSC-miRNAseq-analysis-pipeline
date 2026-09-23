import os
import pytest
import textwrap
from pathlib import Path

FIXTURES = Path(__file__).parent / "fixtures"
SCRIPTS = Path(__file__).parent.parent / "scripts"
REPO_ROOT = Path(__file__).parent.parent


@pytest.fixture
def sample_featurecounts_tsv(tmp_path):
    """Minimal featureCounts output TSV (2 genes, 2 samples)."""
    content = textwrap.dedent("""\
        # Program:featureCounts v2.0.1; Command: ...
        Geneid\tChr\tStart\tEnd\tStrand\tLength\tsample_1.bam\tsample_2.bam
        ENSDARG00000001\tchr1\t1\t1000\t+\t1000\t120\t95
        ENSDARG00000002\tchr1\t2001\t3000\t-\t1000\t340\t410
    """)
    f = tmp_path / "featurecounts.tsv"
    f.write_text(content)
    return str(f)


@pytest.fixture
def sample_gtf(tmp_path):
    """Minimal GTF with gene_biotype and gene_name annotations."""
    content = textwrap.dedent("""\
        chr1\ttest\tgene\t1\t1000\t.\t+\t.\tgene_id "ENSDARG00000001"; gene_name "dre-mir-21"; gene_biotype "miRNA";
        chr1\ttest\tgene\t2001\t3000\t.\t-\t.\tgene_id "ENSDARG00000002"; gene_name "rps11"; gene_biotype "protein_coding";
    """)
    f = tmp_path / "test.gtf"
    f.write_text(content)
    return str(f)


@pytest.fixture
def sample_flagstat(tmp_path):
    """Bowtie2-style flagstat content for one sample."""
    content = textwrap.dedent("""\
        10000 + 0 in total (QC-passed reads + QC-failed reads)
        0 + 0 secondary
        0 + 0 supplementary
        0 + 0 duplicates
        8500 + 0 mapped (85.00% : N/A)
        0 + 0 paired in sequencing
        0 + 0 read1
        0 + 0 read2
        0 + 0 properly paired (N/A : N/A)
        0 + 0 with itself and mate mapped
        0 + 0 singletons (N/A : N/A)
    """)
    f = tmp_path / "sample.flagstat"
    f.write_text(content)
    return str(f)


@pytest.fixture
def sample_bowtie2_log(tmp_path):
    """Bowtie2 alignment log for one sample."""
    content = textwrap.dedent("""\
        10000 reads; of these:
          10000 (100.00%) were unpaired; of these:
            1200 (12.00%) aligned 0 times
            7800 (78.00%) aligned exactly 1 time
            1000 (10.00%) aligned >1 times
        88.00% overall alignment rate
    """)
    f = tmp_path / "sample.bowtie2.mature.log"
    f.write_text(content)
    return str(f)


@pytest.fixture
def sample_bowtie1_log(tmp_path):
    """Bowtie1 alignment log for collapsed reads."""
    content = textwrap.dedent("""\
        # reads processed: 5000
        # reads with at least one alignment: 3800
        # reads that failed to align: 1200
        Reported 3800 alignments
    """)
    f = tmp_path / "sample.bowtie1.hairpin.aln.log"
    f.write_text(content)
    return str(f)


@pytest.fixture
def sample_count_matrix(tmp_path):
    """Small miRNA count matrix CSV for PCA tests (10 miRNAs × 4 samples)."""
    import pandas as pd
    import numpy as np

    rng = np.random.default_rng(42)
    miRNAs = [f"dre-miR-{i}" for i in range(1, 11)]
    samples = ["sample_1", "sample_2", "sample_3", "sample_4"]
    counts = rng.integers(10, 5000, size=(10, 4))
    df = pd.DataFrame(counts, index=miRNAs, columns=samples)
    df.index.name = "miRNA"
    f = tmp_path / "counts.csv"
    df.to_csv(f)
    return str(f)


@pytest.fixture
def sample_isomir_tsv(tmp_path):
    """Minimal miRtop-style isomiR TSV for pivot tests."""
    content = textwrap.dedent("""\
        UID\tRead\tmiRNA\tVariant\tiso_5p\tiso_3p\tiso_add3p\tiso_snp\tsample_1.seqcluster.hairpin.aln.srt\tsample_2.seqcluster.hairpin.aln.srt
        uid1\tTAGCTTATCAGACTGGTGTTGGC\tdre-miR-21\tNA\t0\t0\t0\t0\t120\t95
        uid2\tUAGCUUAUCAGACUGGUGUUGGC\tdre-miR-21\tiso_3p:1\t0\t1\t0\t0\t15\t8
        uid3\tTGAGGTAGTAGGTTGTATAGTT\tdre-let-7a\tNA\t0\t0\t0\t0\t340\t410
    """)
    f = tmp_path / "sample.hairpin.tsv"
    f.write_text(content)
    return str(f)
