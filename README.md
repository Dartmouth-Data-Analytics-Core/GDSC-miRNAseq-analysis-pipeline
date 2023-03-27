# Dartmouth CQB miRNA-seq analysis pipeline
Pipeline for processing and quality control of miRNA-seq data

## Introduction 
This pipeline provides preprocessing and quality control of miRNA sequencing data.  Currently, only miRNA seq runs created with Qiagen chemistry, including UMIs, are supported.  The pipeline has been built and tested using human and mouse data sets. Required software can be installed using Conda with the enrionment file (environment.yml), or specified as paths in the config.yaml file.

## Pipeline summary:
The major steps implmented in the pipeline include: 

- Trimming of adapters and caputuring of UMIs using [*UMI-tools*](https://github.com/CGATOxford/UMI-tools)
- Alignment to mirBase using [*bowtie2*](https://github.com/BenLangmead/bowtie2)
- Alignment of remaining unmapped reads to the whole genome using  [*bowtie2*](https://github.com/BenLangmead/bowtie2)
- Quantification of mirBase miRNAs and other genomic RNAs with [*Samtools*](http://www.htslib.org/) and [*Featurecounts*](http://subread.sourceforge.net/)

All of these tools can be installed in a [conda environment](https://docs.conda.io/en/latest/) or on paths available to a computing server. As input, the pipeline takes raw data in FASTQ format, and produces quantified read counts as well as a quality control report quantifying reads with the correct UMI structure, reads mapping to mirBase, and reads mapping to the genome.

## Implementation
The pipeline uses Snakemake to submit jobs to the scheduler, or spawn processes on a single machine, and requires several variables to be configured by the user when running the pipeline: 
* **sample_tsv** - A TSV file containing sample names and paths to fastq paths.  See example in this repository for formatting.

* **bowtie_index** - Path to mirBase genome reference index  
* **bowtie_genome_index** - Path to bowtie2 genome reference index  

* **annotation_gtf** - Absolute path to genome annotation file (.gtf) of [*Featurecouts*](http://subread.sourceforge.net/) or [*RSEM*](https://deweylab.github.io/RSEM/)
* **featurecounts_strand** - "1" or "2" #1 for first read transcription strand, 2 for second, 0 for unstranded.*  

## Running tests using pre-built environments on Discovery
Clone this repository:
```shell
git clone https://github.com/Dartmouth-Data-Analytics-Core/DAC-miRNAseq-pipeline.git
cd DAC-miRNAseq-pipeline
```
Activate an environment containing Snakemake:
```shell
conda activate /dartfs-hpc/rc/lab/G/GMBSR_bioinfo/misc/sullivan/tools/snakemake/snakemake-7.18
```

## More Command Line Examples
Submit the pipeline to a single machine, allowing usage of 40 cores:
```shell
snakemake --use-conda -s Snakefile -j 40
```

**Contact & questions:** 
Please address questions to *DataAnalyticsCore@groups.dartmouth.edu* or submit an issue in the GitHub repository. 

**This pipeline was created with funds from the COBRE grant **1P20GM130454**. 
If you use the pipeline in your own work, please acknowledge the pipeline by citing the grant number in your manuscript.**
