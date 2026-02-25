# miRNA Pipeline Configuration Parameters

## Overview

This pipeline processes small RNA-seq / miRNA-seq data with configurable options for library layout, spike-ins, UMIs, and reference alignment.

**Alignment strategy:**
- **miRBase alignments** (mature and hairpin miRNAs) are conducted using **Bowtie1**, where **uracils (U) in reads are converted to thymines (T)** prior to alignment.
- **Genome alignment** is conducted using **Bowtie2** against the specified reference genome index.

---

## Configuration Parameters

| Section | Parameter | Allowed Values / Type | Example / Default | Description |
|------|---------|----------------------|------------------|-------------|
| General | `sample_csv` | string (path to CSV) | `sample_fastq_list.csv` | Sample metadata file listing FASTQ files and sample IDs |
| General | `sps` | `"hsa"` or `"dre"` | `hsa` | Species identifier (human or zebrafish) |
| General | `layout` | `"single"` or `"paired"` | `single` | Sequencing layout |
| General | `use_spikeins` | `true` / `false` | `false` | Whether spike-in sequences are used |
| General | `use_umitools` | `true` / `false` | `false` | Whether UMIs are extracted and processed |
| General | `adapter_3prime` | string (DNA sequence) | `AGATCGGAAGAGCACACGTCTGAACTCCAGTCA` | 3′ adapter sequence for cutadapt |
| General | `nextseq_trim` | string (cutadapt flag) or `""` | `""` | NextSeq quality trimming option; empty string disables |

| Spike-in | `sample_with_spikein_finalvolume` | integer | `5` | Final sample volume used for spike-in normalization |
| Spike-in | `bowtie_spikein_index` | string (path) | `/path/to/spikeins_full` | Bowtie index for spike-in reference |
| Spike-in | `spikein_reference_core` | string (path) | `/path/to/spikeins_core.fa` | Core FASTA file for spike-in sequences |

| Software | `umitools_path` | string (path or executable) | `/path/to/umi_tools` | Path to umi-tools binary |
| Software | `samtools_path` | string (path or executable) | `/path/to/samtools` | Path to samtools binary |
| Software | `bowtie1_path` | string (path or executable) | `/path/to/bowtie` | Path to Bowtie1 executable |
| Software | `bowtie2_path` | string (path or executable) | `/path/to/bowtie2` | Path to Bowtie2 executable |
| Software | `multiqc_path` | string (path or executable) | `multiqc` | Path to MultiQC executable |
| Software | `python_mirna_venv` | string (path) | `/path/to/python3` | Python interpreter for miRNA-specific scripts |
| Software | `featurecounts_path` | string (path or executable) | `/path/to/featureCounts` | Path to featureCounts executable |

| Reference (Human) | `bowtie1_mature_index` | string (path) | `/path/to/mature_index` | Bowtie1 index for mature miRNA sequences |
| Reference (Human) | `bowtie1_hairpin_index` | string (path) | `/path/to/hairpins` | Bowtie1 index for hairpin miRNA sequences |
| Reference (Human) | `bowtie2_genome_index` | string (path) | `/path/to/hg38_index` | Bowtie2 genome index |
| Reference (Human) | `annotation_gtf` | string (path) | `gencode.v45.annotation.gtf` | GTF file for genome annotation |
| Reference (Human) | `hairpin_gff` | string (path) | `hsa.gff3` | GFF3 file describing hairpin miRNA loci |
| Reference (Human) | `hairpin_fa` | string (path) | `hsa_mirbase_hairpins.fa` | FASTA file of hairpin miRNA sequences |

| Utilities | `pca_plot_script` | string (path) | `scripts/pca_plotting.py` | Script for PCA visualization |
| Utilities | `featurecounts_strand` | `"0"`, `"1"`, or `"2"` | `0` | Strand specificity for featureCounts (0 = unstranded) |
| Utilities | `featurecounts_annscript` | string (path) | `scripts/add_gene_to_ensg.py` | Script to annotate gene IDs |

---

## Notes

- When `use_umitools` is enabled, UMI extraction and deduplication are applied prior to alignment.
- Spike-in parameters are only used when `use_spikeins: true`.
- Paths may be absolute or resolved via `--use-conda` depending on execution environment.