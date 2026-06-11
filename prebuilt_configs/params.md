# miRNA Pipeline Configuration Parameters

## Overview

This pipeline processes **small RNA-seq / miRNA-seq** data with configurable options for sequencing layout, spike-ins, UMIs, and reference alignment.

### Alignment Strategy

The pipeline performs several alignment and quantification steps:

1. **miRBase alignment (collapsed reads)**  
   Reads are collapsed to miRNA loci using **SeqCluster**, then aligned to a miRNA hairpin index using **Bowtie1**.  
   Uracils in RNA sequences are converted to thymines to ensure compatibility with DNA aligners.

2. **Mature miRNA alignment (uncollapsed reads)**  
   Uncollapsed reads are aligned to a **padded mature miRNA reference** using **Bowtie2**.  
   Unmapped reads from this step are retained for downstream genome alignment.

3. **Genome alignment**  
   Remaining reads are aligned to the full reference genome using **Bowtie2**.

4. **Feature annotation**  
   Genomic alignments are annotated using **featureCounts** with a customized annotation GTF that combines:
   - Ensembl gene annotations
   - Mature tRNA annotations from **GtRNAdb**

---

# Configuration Parameters

## General Parameters

| Parameter | Type / Allowed Values | Example | Description |
|----------|----------------------|--------|-------------|
| `sample_csv` | Path to CSV | `sample_fastq_list.csv` | Sample metadata file listing FASTQ files and sample IDs |
| `sps` | `"hsa"` or `"dre"` or `"mmu"` | `hsa` | Species identifier (human or zebrafish or mouse) |
| `layout` | `"single"` or `"paired"` | `single` | Sequencing layout |
| `use_spikeins` | `true` / `false` | `false` | Whether spike-in sequences are used |
| `use_umitools` | `true` / `false` | `false` | Whether UMIs are extracted and processed |
| `adapter_3prime` | DNA sequence | `AGATCGGAAGAGCACACGTCTGAACTCCAGTCA` | 3′ adapter sequence used by cutadapt |
| `nextseq_trim` | cutadapt flag or empty string | `""` | NextSeq-specific quality trimming option |

---

## Spike-in Parameters

Used only when `use_spikeins: true`.

| Parameter | Type | Example | Description |
|----------|------|--------|-------------|
| `sample_with_spikein_finalvolume` | Integer | `8` | Final sample volume used for spike-in normalization |
| `bowtie_spikein_index` | Path | `/path/to/spikeins_full` | Bowtie index for spike-in reference |
| `spikein_reference_core` | Path | `/path/to/spikeins_core.fa` | FASTA file containing core spike-in sequences |

---

## Software Paths

These parameters define paths to required executables.

| Parameter | Type | Example | Description |
|----------|------|--------|-------------|
| `umitools_path` | Path or executable | `/path/to/umi_tools` | Path to umi-tools |
| `samtools_path` | Path or executable | `/path/to/samtools` | Path to samtools |
| `bowtie1_path` | Path or executable | `/path/to/bowtie` | Path to Bowtie1 |
| `bowtie2_path` | Path or executable | `/path/to/bowtie2` | Path to Bowtie2 |
| `multiqc_path` | Path or executable | `multiqc` | Path to MultiQC |
| `python_mirna_venv` | Path | `/path/to/python3` | Python interpreter used for miRNA scripts |
| `featurecounts_path` | Path or executable | `/path/to/featureCounts` | Path to featureCounts |

---

## Reference Files (Human)

| Parameter | Type | Example | Description |
|----------|------|--------|-------------|
| `bowtie1_mature_index` | Path | `/path/to/mature_index` | Bowtie1 index for mature miRNA sequences |
| `bowtie1_hairpin_index` | Path | `/path/to/hairpins` | Bowtie1 index for miRNA hairpin sequences |
| `bowtie2_genome_index` | Path | `/path/to/hg38_index` | Bowtie2 genome index |
| `annotation_gtf` | Path | `hg38_smRNAs_super_sorted.gtf` | Custom gtf, genome + tRNA information |
| `hairpin_gff` | Path | `hsa.gff3` | GFF3 file describing miRNA hairpin loci |
| `hairpin_fa` | Path | `hsa_mirbase_hairpins.fa` | FASTA file of miRNA hairpin sequences |

---

## Utility Parameters

| Parameter | Type | Example | Description |
|----------|------|--------|-------------|
| `featurecounts_strand` | `"0"`, `"1"`, `"2"` | `0` | Strand specificity for featureCounts (0 = unstranded) |

---

# Notes

- When `use_umitools` is enabled, UMI extraction and deduplication are applied prior to alignment.
- Spike-in parameters are ignored unless `use_spikeins` is enabled.
- Software paths can point to system executables or Conda-managed environments.
- Reference paths should point to pre-built indexes compatible with the specified aligners.
