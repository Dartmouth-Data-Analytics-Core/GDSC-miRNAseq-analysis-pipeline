#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# GDSC miRNA Pipeline v3
#
# Pipeline for the quantification of miRNAs, isomiRs, and other small RNAs
#
# TO DO
# - Add script to collapse isomirs down to their family to provide a decent proxy for mature counts
#
# TO-DO
#--------
# - From clover-seq gtf file, remove miRNA rows
# - Find piRNA gff3
#   - Convert gff3 to gtf
#   - Find circRNA gtf
#   - Append these to human genome gtf (try and remove multiple records?)
#   - Align to mirBase mature, unaligned get mapped to human genome
#   - Annotate genome hits with feature counts with "super gtf"
# 
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
import pandas as pd
import pprint

#----- set config file if not defined in job script
configfile: "config.yaml"
USE_SPIKEINS = config.get("use_spikeins", False)
USE_UMITOOLS = config.get("use_umitools", False)

#----- read in sample data
samples_df = pd.read_csv(config["sample_csv"]).set_index("sample_id", drop=False)
sample_list = list(samples_df['sample_id'])

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# RULE ALL INPUTS
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#----- Include additional rules
include: "additional_rules/genome_alignment/genome_aln.smk"
include: "additional_rules/QC/qc.smk"
if USE_SPIKEINS:
    include: "additional_rules/spike_ins/spikein.smk"
if USE_UMITOOLS:
    include: "additional_rules/umitools/umi_extract.smk"

#----- Build main pipeline input list
all_inputs = []

#----- Trimming
all_inputs += expand("trimming/{sample}.R1.trim.fastq.gz", sample=sample_list)

#----- QC
all_inputs += ["fastQC/fastqc_multiqc_config.yaml"]
all_inputs += expand("fastQC/{sample}.R1.trim_fastqc.html", sample=sample_list)
all_inputs += expand("fastQC/{sample}.R1.trim_fastqc.zip", sample=sample_list)

#----- UMI deduplication outputs (optional)
if USE_UMITOOLS:
    all_inputs += expand("umi_reads/{sample}.umi.fastq.gz", sample=sample_list)

#----- Seqcluster collapsed reads
all_inputs += expand("collapsed/{sample}.seqcluster.fastq.gz", sample=sample_list)
all_inputs += expand("collapsed/{sample}.seqcluster.hairpin.aln.srt.bam", sample=sample_list)

#----- miRBase alignment and metrics (always included)
all_inputs += expand("mirbase_alignment/{sample}.mature.srt.bam", sample=sample_list)
if USE_UMITOOLS:
    all_inputs += expand("mirbase_alignment/{sample}.mature.srt.dedup.bam", sample=sample_list)
all_inputs += expand("mirbase_alignment/{sample}.mature.unalign.fastq", sample=sample_list)
all_inputs += expand("mirbase_alignment/{sample}.mature.srt.bam.idxstats", sample=sample_list)
all_inputs += expand("mirbase_alignment/{sample}.mature.srt.bam.flagstat", sample=sample_list)
all_inputs += ["miRNA_Quant/raw_merged_canonical_and_all_isomirs.csv"]

#----- mirtop outputs
all_inputs += expand("mirtop/{sample}.hairpin.gff", sample=sample_list)
all_inputs += expand("mirtop/temp/{sample}.hairpin_long.csv", sample=sample_list)
all_inputs += ["mirtop/mirtop_stats.log"]

#----- Genome alignment and featureCounts (always included)
all_inputs += expand("genome_alignment/{sample}.genome.srt.filt.bam", sample=sample_list)
if USE_UMITOOLS:
    all_inputs += expand("genome_alignment/{sample}.genome.srt.filt.dedup.bam", sample=sample_list)
all_inputs += [
    "genome_counts/featurecounts.tsv",
    "genome_counts/featurecounts.readcounts.tsv",
    "genome_counts/featurecounts.readcounts.biotype.tsv",
    "metrics/mirna_genome_alignment_metrics.tsv",
    "metrics/mirna_genome_alignment_metrics.xlsx"]

#----- PCA
all_inputs += [
    "plots/PCA_top_PC1_vs_PC2.png",
    "plots/PCA_top_PCA_variance_bar.png"]

#----- Spike-in outputs (optional)
if USE_SPIKEINS:
    all_inputs += expand("spikein_alignment/{sample}.unmapped.bowtie.fastq.gz", sample=sample_list)
    all_inputs += [
        "spikein_counts/spikein.readcounts.tsv",
        "spikein_metrics/normalized_scalefactor_canon_and_isomir_counts.tsv"]

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# PIPELINE
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#----- Main pipieline execution
rule all:
    input:
        all_inputs
    output:
        "multiqc_report.html",
        "fastQC/trimmed_fastqc_report.html"       
    conda:
        "env_config/multiqc.yaml",
    resources: cpus="10", maxtime="2:00:00", mem_mb="60gb",
    params:
        multiqc=config["multiqc_path"],
        use_umi = USE_UMITOOLS,
    shell: """

        #----- Run multiqc on fastq results
        multiqc \
            fastQC \
            -c fastQC/fastqc_multiqc_config.yaml
        mv multiqc_report.html fastQC/trimmed_fastq_report.html &&
        mv multiqc_data fastQC/multiqc_data &&

        #----- Run multiqc
        if [ "{params.use_umi}" = "true" ]; then
            {params.multiqc} -v -c multiqc_config.yaml \
                alignment_logs/mirbase_mature_padded \
                alignment_logs/genome_alignment \
                mirbase_alignment \
                genome_counts \
                umi_reads
        else
            {params.multiqc} -v -c multiqc_config.yaml \
                alignment_logs/mirbase_mature_padded \
                alignment_logs/genome_alignment \
                mirbase_alignment \
                genome_counts \
                mirtop
        fi

        #----- Clean
        if [ -d "mirtop/log" ];
            rm -r mirtop/log
        fi

"""

#----- Rule to execute trimming
rule trimming:
    """
    Read trimming
    """
    output: 
        "trimming/{sample}.R1.trim.fastq.gz",
    params:
        sample = lambda wildcards:  wildcards.sample,
        fastq_file_1 = lambda wildcards: samples_df.loc[wildcards.sample, "fastq_1"],
        adapter_3prime = config["adapter_3prime"],
        nextseq_trim = config["nextseq_trim"],
    conda:
        "env_config/cutadapt.yaml",
    resources: 
        cpus="10", 
        maxtime="2:00:00", 
        mem_mb="60gb",
    message: "Trimming {wildcards.sample} reads with cutadapt."
    shell: """

        #----- Make log directory
        mkdir -p trimming/logs

        #----- Run cutadapt
        cutadapt \
            -o trimming/{params.sample}.R1.trim.fastq.gz \
            {params.fastq_file_1} \
            -m 1 \
            {params.nextseq_trim} \
            -j {resources.cpus} \
            -q 30 \
            --max-n 0.8 \
            -a {params.adapter_3prime} \
            --trim-n > trimming/logs/{params.sample}.cutadapt.log
    """

#----- Function selecting input FASTQ file for alignment
def get_alignment_input(wildcards):
    if USE_UMITOOLS:
        return f"umi_reads/{wildcards.sample}.umi.fastq.gz"
    return f"trimming/{wildcards.sample}.R1.trim.fastq.gz"

#----- Rule to collapse reads with seqcluster
rule seqcluster:
    """
    Collapse reads
    """
    input: get_alignment_input
    output: 
        collapsed = "collapsed/{sample}.seqcluster.fastq.gz"
    params:
        sample = lambda wildcards: wildcards.sample,
        seqcluster_path = config["seqcluster_path"]
    conda: "env_config/seqcluster.yaml"
    threads: 8
    resources:
        maxtime="2:00:00",
        mem_mb="60gb"
    message: "Collapsing {wildcards.sample} reads with Seqcluster."
    log: "collapsed/logs/{sample}.seqcluster.log"
    shell: """
    
        #----- Make log subdirectory
        mkdir -p collapsed/logs

        #----- Run seqcluster to collapse reads
        seqcluster \
            collapse \
            -m 1 \
            --min_size 15 \
            -f {input} \
            -o collapsed > {log} 
        
        #----- Clean
        mv collapsed/{params.sample}.R1.trim_trimmed.fastq collapsed/{params.sample}.seqcluster.fastq
        gzip collapsed/{params.sample}.seqcluster.fastq
    
    """

#----- Rule to align collapsed reads to hairpins
rule collapsed_hairpin_aln:
    """
    Align collapsed reads to hairpin sequences
    """
    input:
        seqcluster_fastq = "collapsed/{sample}.seqcluster.fastq.gz"
    output:
        collapsed_aln = "collapsed/{sample}.seqcluster.hairpin.aln.srt.bam"
    params:
        sample = lambda wildcards:  wildcards.sample,
        bowtie1_path = config["bowtie1_path"],
        hairpin_index = config["bowtie1_hairpin_index"],
        samtools_path = config["samtools_path"]
    threads: 8
    resources:
        maxtime="2:00:00",
        mem_mb="60gb"
    message: "Aligning {wildcards.sample} collapsed reads to hairpin index with Bowtie1."
    log: "alignment_logs/seqcluster/{sample}.bowtie1.hairpin.aln.log"
    shell: """
    
        #----- Make logs subdirectory
        mkdir -p alignment_logs/seqcluster

        #----- Align collapsed reads to hairpins
        {params.bowtie1_path} \
            --threads {threads} \
            --sam \
            -x {params.hairpin_index} \
            -q \
            --un collapsed/{params.sample}.unmapped.fastq \
            -t \
            -k 50 \
            --best \
            --strata \
            -e 99999 \
            --chunkmbs 2048 \
            {input.seqcluster_fastq} \
            2>| >(tee {log} >&2) \
            | {params.samtools_path} view -@ 24 -bS - \
            | {params.samtools_path} sort -@ 24 -o {output.collapsed_aln}

        {params.samtools_path} index {output.collapsed_aln}
    
    """

#----- Rule to run miRtop
rule miRtop:
    """
    Running mirtop for isomiR calculation
    """
    input:
        collapsed_aln = "collapsed/{sample}.seqcluster.hairpin.aln.srt.bam"
    output:
        mirtop_gff = "mirtop/{sample}.hairpin.gff",
        hairpin_tsv = "mirtop/{sample}.hairpin.tsv"
    conda: "env_config/mirtop.yaml"
    params:
        sample = lambda wildcards:  wildcards.sample,
        hairpin_fa = config["hairpin_fa"],
        hairpin_gff = config["hairpin_gff"],
        sps = config["sps"]
    resources:
        cpus="10", 
        maxtime="2:00:00", 
        mem_mb="60gb",
    message: "Getting {wildcards.sample} isomiRs with miRtop."
    log: 
        gffLog = "mirtop/logs/{sample}.mirtop.gff.log",
        countLog = "mirtop/logs/{sample}.mirtop.counts.log",
        statLog = "mirtop/logs/{sample}.mirtop.stats.log"
    shell: """

        #----- Make logs subdirectory
        mkdir -p mirtop/logs

        #----- Run miRtop GFF
        mirtop gff \
            --add-extra \
            --sps {params.sps} \
            --hairpin {params.hairpin_fa} \
            --gtf {params.hairpin_gff} \
            -o mirtop \
            {input.collapsed_aln} > {log.gffLog} 2>&1 &&
        mv mirtop/{params.sample}.seqcluster.hairpin.aln.srt.gff mirtop/{params.sample}.hairpin.gff
        
        #----- Run miRtop counts
        mirtop counts \
            -o mirtop \
            --hairpin {params.hairpin_fa} \
            --gff {output.mirtop_gff} \
            --gtf {params.hairpin_gff} > {log.countLog} 2>&1

"""

#----- Rule to run mirtop stats
rule mirtop_stats:
    """
    Run mirtop stats
    """
    input:
        expand("mirtop/{sample}.hairpin.gff", sample=sample_list)
    output:
        "mirtop/mirtop_stats.log"
    conda: "env_config/mirtop.yaml"
    resources:
        cpus="10", 
        maxtime="2:00:00", 
        mem_mb="60gb",
    message: "Getting mirtop stats"
    shell: """

        #----- Run mirtop stats
        mirtop stats \
            {input} \
            -o mirtop
    
    """

#----- Rule to pivot isomiRs longer
rule pivot_isomirs_longer:
    """
    Pivot isomiRs longer
    """
    input:
        isomiR_counts = "mirtop/{sample}.hairpin.tsv"
    output:
        isomiR_long = "mirtop/temp/{sample}.hairpin_long.csv"
    conda: "env_config/r_env.yaml"
    params:
        sample = lambda wildcards:  wildcards.sample,
    resources:
        cpus="10", 
        maxtime="2:00:00", 
        mem_mb="60gb",
    message: "Pivottings {wildcards.sample} isomiR data to long format."
    shell: """

        #----- Pivot longer
        Rscript scripts/pivot_longer.R \
            {input} \
            {output}
    
    """

#----- Rule to collate master isomiR table, remove canonical miRNA counts from isomiR table
rule collate_isomir_table:
    """
    Combine isomiR results into master table
    """
    input:
        expand("mirtop/temp/{sample}.hairpin_long.csv", sample=sample_list)
    output:
        "miRNA_Quant/raw_merged_canonical_and_all_isomirs.csv"
    conda: "env_config/r_env.yaml"
    resources:
        cpus="10", 
        maxtime="2:00:00", 
        mem_mb="60gb",
    message: "Collating isomir counts across samples."
    shell: """
    
    #----- Collate all counts files
    awk 'NR == 1 || FNR > 1' \
        {input} > "mirtop/temp/master_counts_long.csv"

    #----- Create master counts file with formatting
    Rscript scripts/pivot_wider.R \
        mirtop/temp/master_counts_long.csv \
        miRNA_Quant/
    
    
    """

#----- Rule to run PCA
rule pca_plots:
    """
    Run PCA on the isomir data
    """
    input: 
        "miRNA_Quant/raw_merged_canonical_and_all_isomirs.csv",
    output:
        "plots/PCA_top_PC1_vs_PC2.png",
        "plots/PCA_top_PCA_variance_bar.png",
    conda:
        "env_config/pcaplot.yaml",
    params:
        pca_plot_script = config['pca_plot_script'],   
    resources: cpus="1", maxtime="1:00:00", mem_mb=2000,
    message: "Running PCA"
    shell: """

        #----- Create directory
        mkdir -p plots

        #----- Run PCA script
        python scripts/pca_plotting.py \
            {input} \
            plots
    """





