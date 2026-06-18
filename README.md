# Dartmouth GDSC miRNA-seq analysis pipeline
<img src="img/cqb_logo.jpg" alt="CQB Logo" width="200" align="right"/>  

![Version](https://img.shields.io/badge/version-2.0-blue)

The GDSC miRNA-seq pipeline provides the preprocessing of microRNA-Seq (miRNA) data with robust quality control and data visualization implemented through [Snakemake](https://snakemake.readthedocs.io/en/stable/) for use on the [Dartmouth Discovery HPC](https://rc.dartmouth.edu/discoveryhpc/). This pipeline supports the quantification of isomiRs (miRNA sequence variants) and canonical mature miRNAs via [miRBase](https://www.mirbase.org). Currently, this pipeline is highly configurable and compatible with human (hg38), mouse (mm10), and zebrafish (grcz11) across non-UMI containing Qiagen libraries and UMI-containing NEB Small RNA chemistries. Robust testing was conducted on human and zebrafish data. Required software can be installed using Conda.

## Pipeline summary:

**The major steps implmented in the pipeline include:**

- Trimming of adapters using [*Cutadapt*](https://github.com/marcelm/cutadapt). 
- Collapsing of trimmed reads to miRNA loci using [*seqcluster*](https://github.com/lpantano/seqcluster). 
- Alignment of collapsed reads to mirBase hairpin sequences using [*bowtie1*](https://github.com/BenLangmead/bowtie). 
- IsomiR quantification (canonical miRNAs + isomiRs) using [*miRTop*](https://github.com/miRTop/mirtop). 
- Alignment of trimmed un-collapsed reads to padded mature miRbase sequence using [*bowtie2*](https://github.com/BenLangmead/bowtie2). 
- Alignment of unaligned reads to the full genome using [*bowtie2*](https://github.com/BenLangmead/bowtie2). 
- Quantification of genomic sequences and other small RNAs (tRNAs, snoRNAs, Mt-rRNA, Mt-tRNAs, etc.) with [*Samtools*](http://www.htslib.org/) and [*Featurecounts*](http://subread.sourceforge.net/).

**Optional pipeline steps include:**

- Capturing of unique molecular identifiers (UMIs) for deduplication using [*UMI-tools*](https://github.com/CGATOxford/UMI-tools). 
- Spike-in data support and normalization using [*BBMap*](https://github.com/BioInfoTools/BBMap) and [*bowtie2*](https://github.com/BenLangmead/bowtie2).

All of these tools can be installed in a [conda environment](https://docs.conda.io/en/latest/) or on paths available to a computing server. As input, the pipeline takes raw data in FASTQ format, and produces quantified read counts as well as a quality control report quantifying reads with the correct UMI structure, isomiR counts, canonical miRNA counts, and a merged isomiR + canonical counts (grouped by miRNA), reads mapping to the genome, and associated genome counts with other small RNA biotypes.

## Implementation

The pipeline uses Snakemake to submit jobs to the scheduler, or spawn processes on a single machine, and requires several variables to be configured by the user when running the pipeline. For ease, prebuilt configuration files have been built for human and zebrafish.

> [!IMPORTANT]
> For a detailed description of available parameters, see [`prebuilt_configs/params.md`](prebuilt_configs/params.md)


* **sample_csv** - A CSV file containing sample names and paths to fastq paths.  See example in this repository for formatting.

* **bowtie1_hairpin_index** - Path to bowtie1 mirBase hairpin index, built and formatted for compatibility with miRTop.  
* **hairpin_gff** - Path to GFF3 annotation for hairpin sequences, used for isomiR analysis with miRTop.  
* **hairpin_fa** - Path to FASTA file of miRNA hairpin sequences used for isomiR analysis with miRTop.

* **padded_mature_index** - Path to bowtie2 mirBase mature index. 
* **bowtie2_genome_index** - Path to bowtie2 genome index.  
* **annotation_gtf** - Path to modified organism gtf for featurecounts
* **featurecounts_strand** - "1" or "2" #1 for first read transcription strand, 2 for second, 0 for unstranded.*    
* **adapter_3prime** - the 3'-end adapter that was used in sequencing. It is important this is specified correctly to ensure fidelity of end-to-end alignment, which will not soft clip sequencing adapters.   
* **nextseq_trim** - cutadapt flag and value for NextSeq quality trimming; set to "" to disable. 

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

Modify the `job.script.sh` to point to a valid organism to dynamically set a prebuilt config:

```shell
#----- Specify Config (one of "human", "mouse", or "zebrafish", case-sensitive and needs to be in quotes.)
CONFIG="human"
```

Submit the `job.script.sh`

```shell
sbatch job.script.sh
```


**Contact & questions:** 
Please address questions to *DataAnalyticsCore@groups.dartmouth.edu* or submit an issue in the GitHub repository. 

**This pipeline was created with funds from the COBRE grant **1P20GM130454**. 
If you use the pipeline in your own work, please acknowledge the pipeline by citing the grant number in your manuscript.**
