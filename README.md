# Dartmouth GDSC miRNA-Seq Pipeline
<img src="img/cqb_logo.jpg" alt="CQB Logo" width="200" align="right"/>

![Version](https://img.shields.io/badge/version-2.0-blue)

The GDSC miRNA-Seq pipeline provides preprocessing and quantification of microRNA sequencing data with robust quality control and data visualization, implemented through [Snakemake](https://snakemake.readthedocs.io/en/stable/) for use on the [Dartmouth Discovery HPC](https://rc.dartmouth.edu/discoveryhpc/). The pipeline supports quantification of both canonical mature miRNAs and isomiRs (miRNA sequence variants) via [miRBase](https://www.mirbase.org), and is compatible with human (hg38), mouse (mm10), and zebrafish (GRCz11) across non-UMI Qiagen libraries and UMI-containing NEB Small RNA chemistries.

## Documentation

- [Summary](#summary)
- [Installation](#installation)
- [Configuration](#configuration)
- [Optional Features](#optional-features)
- [Reference Integrity](#reference-integrity)
- [Parameters](prebuilt_configs/params.md)
- [Understanding the Outputs](understanding_outputs.md)
- [Contact](#contact)

## Summary

The pipeline supports the use of Conda environments for all software dependencies as well as singularity containers. Environment files are located in [`env_config/`](env_config/). Singularity images are hosted by the GDSC. 

To run this pipeline:
1. Populate [`sample_fastq_list.csv`](sample_fastq_list.csv) with your sample information
2. Set the `CONFIG` variable in [`job.script.sh`](job.script.sh) to your organism
3. Adjust any parameters in [`prebuilt_configs/`](prebuilt_configs/) if needed (see [Parameters](prebuilt_configs/params.md))
4. Submit [`job.script.sh`](job.script.sh) to the SLURM scheduler

Currently the pipeline performs the following:

- Adapter trimming with [Cutadapt](https://cutadapt.readthedocs.io/en/stable/)
- Collapsing of trimmed reads to miRNA loci using [seqcluster](https://github.com/lpantano/seqcluster)
- Alignment of collapsed reads to miRBase hairpin sequences using [Bowtie1](https://github.com/BenLangmead/bowtie)
- IsomiR quantification (canonical miRNAs + isomiRs) using [miRTop](https://github.com/miRTop/mirtop)
- Alignment of trimmed reads to padded mature miRBase sequences using [Bowtie2](https://github.com/BenLangmead/bowtie2)
- Alignment of unaligned reads to the full genome using [Bowtie2](https://github.com/BenLangmead/bowtie2)
- Quantification of genomic features and small RNA biotypes (tRNAs, snoRNAs, Mt-rRNA, Mt-tRNA, etc.) using [Samtools](http://www.htslib.org/) and [featureCounts](http://subread.sourceforge.net/)
- Quality control and summary reporting using [FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/) and [MultiQC](https://multiqc.info/)
- PCA and variance plots from isomiR count matrices using a custom Python script
- **(optional)** UMI extraction and deduplication using [UMI-tools](https://github.com/CGATOxford/UMI-tools)
- **(optional)** Spike-in alignment and normalization using [BBMap](https://github.com/BioInfoTools/BBMap) and [Bowtie](https://github.com/BenLangmead/bowtie)

## Installation

Clone the repository:

```shell
git clone https://github.com/Dartmouth-Data-Analytics-Core/GDSC-miRNAseq-analysis-pipeline
cd GDSC-miRNAseq-analysis-pipeline
```

## Configuration

**1. Sample sheet**

Populate [`sample_fastq_list.csv`](sample_fastq_list.csv) with your sample information. This is a comma-separated file with the following columns:

| Column | Description |
|--------|-------------|
| `sample_id` | Short sample identifier used to name all output files |
| `fastq_1` | Path to the R1 FASTQ file |

**2. Job submission script**

>[!IMPORTANT]
> You must set the `CONFIG` variable in [`job.script.sh`](job.script.sh) before submitting. Accepted values are `human`, `mouse`, or `zebrafish` (case-sensitive).

```shell
CONFIG="human"
```

Setting `CONFIG` automatically selects the correct prebuilt config file (`prebuilt_configs/${CONFIG}_config.yaml`). No other path changes are required when using a prebuilt config.

**3. Pipeline parameters**

Each organism has a prebuilt config in [`prebuilt_configs/`](prebuilt_configs/). These files contain all tunable settings for trimming, alignment, and optional analyses. See [Parameters](prebuilt_configs/params.md) for a full description of every parameter. For an explanation of all output files, see [Understanding the Outputs](understanding_outputs.md).

**4. Reference integrity**

All reference files (miRBase indices, genome indices, annotation GTF/GFF) were built from [miRBase v22](https://www.mirbase.org/blog/2018/11/mirbase-22-released/). These references are intended to be static, but the pipeline includes an MD5 checksum check at startup to guard against unintended reference drift (e.g. accidental overwrites or storage corruption). Pre-computed checksums for each organism are stored in [`prebuilt_configs/ref_md5s/`](prebuilt_configs/ref_md5s/). **DO NOT EDIT THESE FILES BY HAND.**

**5. Submitting the job**

```shell
sbatch job.script.sh
```

## Optional Features

### UMI-tools deduplication

>[!IMPORTANT]
> UMI support is designed for NEB Small RNA libraries containing a 12-base UMI. Set `use_umitools: true` in your config to enable UMI extraction and deduplication.

```yaml
use_umitools: true
```

When enabled, UMIs are extracted from raw reads prior to trimming using [UMI-tools](https://github.com/CGATOxford/UMI-tools), and PCR duplicates are removed after alignment. Deduplicated BAMs are used for all downstream quantification steps.

### Spike-in normalization

>[!IMPORTANT]
> Spike-in support is designed for NextFLEX libraries. Set `use_spikeins: true` and configure `sample_with_spikein_finalvolume` in your config before running.

```yaml
use_spikeins: true
sample_with_spikein_finalvolume: 8
```

When enabled, reads are first aligned to spike-in sequences using [BBMap](https://github.com/BioInfoTools/BBMap) and [Bowtie](https://github.com/BenLangmead/bowtie). Unmapped reads proceed through the standard pipeline. Spike-in counts are used to compute normalization scale factors applied to the final isomiR count matrices.

## Contact

**Contact and questions:** Please address questions to *DataAnalyticsCore@groups.dartmouth.edu* or submit an issue in the GitHub repository.

**This pipeline was created with funds from the COBRE grant 1P20GM130454. If you use the pipeline in your own work, please acknowledge the pipeline by citing the grant number in your manuscript.**
