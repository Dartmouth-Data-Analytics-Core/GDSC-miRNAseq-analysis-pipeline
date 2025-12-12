#####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# setup environment
#####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
import pandas as pd

# set config file
configfile: "config.yaml"
USE_SPIKEINS = config.get("use_spikeins", False)
USE_UMITOOLS = config.get("use_umitools", False)

# read in sample data
samples_df = pd.read_table(config["sample_tsv"]).set_index("sample_id", drop=False)
sample_list = list(samples_df['sample_id'])

#####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# define rules
#####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

if USE_SPIKEINS:
    include: "additional_rules/spike_ins/spikein.smk"
if USE_UMITOOLS:
    include: "additional_rules/umitools/umi_extract.smk"

rule all:
    input:
        expand("trimming/{sample}.R1.trim.fastq.gz", sample=sample_list),
        expand("trimming/{sample}.cutadapt.report", sample=sample_list),    
        expand("umi_reads/{sample}.umi.fastq.gz", sample=sample_list) if USE_UMITOOLS else [],
        expand("mirbase_alignment/{sample}.srt.bam", sample=sample_list),
        expand("mirbase_alignment/{sample}.srt.dedup.bam", sample=sample_list) if USE_UMITOOLS else [],
        expand("mirbase_alignment/{sample}.unalign.fastq", sample=sample_list),
        expand("genome_alignment/{sample}.srt.bam", sample=sample_list),
        expand("genome_alignment/{sample}.srt.filt.bam", sample=sample_list),
        expand("genome_alignment/{sample}.srt.filt.dedup.bam", sample=sample_list) if USE_UMITOOLS else [],
        "metrics/mirna_genome_alignment_metrics.tsv",
        "genome_counts/featurecounts.readcounts.ann.tsv",
        "genome_counts/featurecounts.readcounts_tpm.tsv",
        "genome_counts/featurecounts.readcounts_tpm.ann.tsv",
        "mirbase_counts/mirbase.readcounts.tsv",
        "mirbase_counts/mirbase.readcounts_tpm.tsv",
        expand("spikein_alignment/{sample}.unmapped.bowtie.fastq.gz", sample=sample_list) if USE_SPIKEINS else [],
        "spikein_counts/spikein.readcounts.tsv" if USE_SPIKEINS else [],
        expand("mirbase_alignment/{sample}.srt.bam.idxstats", sample=sample_list),
        expand("mirbase_alignment/{sample}.srt.bam.flagstat", sample=sample_list),
        "spikein_metrics/normalized_scalefactor_mirbase_counts.tsv" if USE_SPIKEINS else [],
        "plots/PCA_1_vs_2.png"
                
    conda:
        "env_config/multiqc.yaml",
    resources: cpus="10", maxtime="2:00:00", mem_mb="60gb",

    params:
        multiqc=config["multiqc_path"],
        use_umi = USE_UMITOOLS,
    output:
        "multiqc_report.html"

    shell: """
        if [ "{params.use_umi}" = "true" ]; then
            {params.multiqc} -v -c multiqc_config.yaml genome_alignment mirbase_alignment genome_counts mirbase_counts umi_reads
        else
            {params.multiqc} -v -c multiqc_config.yaml genome_alignment mirbase_alignment genome_counts mirbase_counts
        fi
"""

rule trimming:
    output: 
        "trimming/{sample}.R1.trim.fastq.gz",
        "trimming/{sample}.cutadapt.report"
    params:
        sample = lambda wildcards:  wildcards.sample,
        fastq_file_1 = lambda wildcards: samples_df.loc[wildcards.sample, "fastq_1"],
    conda:
        "env_config/cutadapt.yaml",
    resources: cpus="10", maxtime="2:00:00", mem_mb="60gb",
    shell: """
        cutadapt \
            -o trimming/{params.sample}.R1.trim.fastq.gz \
            {params.fastq_file_1} \
            -m 1 \
            --nextseq-trim=30 \
            -j {resources.cpus} \
            -q 30 \
            --max-n 0.8 \
            --trim-n > trimming/{params.sample}.cutadapt.report
    """

# define function selecting input FASTQ file for alignment
def get_alignment_input(wildcards):
    if USE_UMITOOLS:
        return f"umi_reads/{wildcards.sample}.umi.fastq.gz"
    return f"trimming/{wildcards.sample}.R1.trim.fastq.gz"

rule mirbase_alignment:
    input:
        get_alignment_input
    output:
        "mirbase_alignment/{sample}.srt.bam",
        "mirbase_alignment/{sample}.unalign.fastq",
    params:
        sample = lambda wildcards:  wildcards.sample,
        bowtie_path = config["bowtie_path"],
        bowtie_index = config["bowtie_index"],
        samtools_path = config["samtools_path"],

    resources: cpus="10", maxtime="2:00:00", mem_mb="60gb",

    shell: """
          {params.bowtie_path} \
            -x {params.bowtie_index} \
            -U {input} -p 12 \
            --norc \
            -D 20 -R 3 -N 1 -L 12 -i S,1,0.50 \
            --un mirbase_alignment/{params.sample}.unalign.fastq \
            -S mirbase_alignment/{params.sample}.aln.sam 2>mirbase_alignment/{params.sample}_mirbase.log.txt

        # subset reads for aligned length > 16 & < 28bp & any reads with gaps (XO/XG tags)
        {params.samtools_path} view -h mirbase_alignment/{params.sample}.aln.sam | \
            awk 'BEGIN {{OFS="\t"}} $1 ~ /^@/ || ((length($10) > 16 && length($10) <= 28) && ($0 !~ /XG:i:[^0]/ && $0 !~ /XO:i:[^0]/)) {{print $0}}' | \
            samtools view -Sb - > mirbase_alignment/{params.sample}.bam
        # filter for any reads with MAPQ <=1
        {params.samtools_path} view -h -q 2 mirbase_alignment/{params.sample}.bam > mirbase_alignment/{params.sample}.sub.bam
        # filter for any reads with > 2 mismatches 
        {params.samtools_path} view -h mirbase_alignment/{params.sample}.sub.bam | \
            awk 'BEGIN {{OFS="\t"}} /^@/ || ($0 ~ /NM:i:[0-2]($|\t)/)' | \
            samtools view -b > mirbase_alignment/{params.sample}.sub2.bam
        
        # sort and index BAM file 
        {params.samtools_path} sort -@ 4 mirbase_alignment/{params.sample}.sub2.bam > mirbase_alignment/{params.sample}.srt.bam
        {params.samtools_path} index mirbase_alignment/{params.sample}.srt.bam

        # remove intermediate bam files 
        rm -rf mirbase_alignment/{params.sample}.bam
        rm -rf mirbase_alignment/{params.sample}.sub.bam
        rm -rf mirbase_alignment/{params.sample}.sub2.bam 
"""

# define function selecting input BAM file for mirbase_stats
def get_mirbase_stats_input(wildcards):
    if USE_UMITOOLS:
        return f"mirbase_alignment/{wildcards.sample}.srt.dedup.bam"
    return f"mirbase_alignment/{wildcards.sample}.srt.bam"

rule mirbase_stats:
    input: 
        get_mirbase_stats_input,
    output:
        "mirbase_alignment/{sample}.srt.bam.idxstats",
        "mirbase_alignment/{sample}.srt.bam.flagstat"
    params:
        sample = lambda wildcards:  wildcards.sample,
        bowtie_path = config["bowtie_path"],
        bowtie_index = config["bowtie_index"],
        samtools_path = config["samtools_path"],

    resources: cpus="10", maxtime="2:00:00", mem_mb="60gb",

    shell: """
    {params.samtools_path} idxstats {input} > mirbase_alignment/{params.sample}.srt.bam.idxstats
    {params.samtools_path} flagstat {input} > mirbase_alignment/{params.sample}.srt.bam.flagstat

"""


rule mirbase_count:
    input:
        expand("mirbase_alignment/{sample}.srt.bam.idxstats", sample=sample_list),

    output:
        "mirbase_counts/mirbase.readcounts.tsv",
        "mirbase_counts/mirbase.readcounts_tpm.tsv",
    params:

    resources: cpus="10", maxtime="2:00:00", mem_mb="60gb",

    shell: """
    echo -ne mirbase_ID"\t"Length"\t" > mirbase_counts/mirbase.readcounts.tsv
    echo {input} | tr " " "\t"| sed s/"mirbase_alignment\/"//g| sed s/".srt.bam.idxstats"//g >> mirbase_counts/mirbase.readcounts.tsv
    paste {input}| awk -f scripts/mirbase_counts.awk >> mirbase_counts/mirbase.readcounts.tsv

    # run TPM normalization 
    python scripts/mirbase-readcnt_to_tpm.py mirbase_counts/mirbase.readcounts.tsv
"""


rule genome_alignment:
    input: 
        "mirbase_alignment/{sample}.unalign.fastq",
    output:
        "genome_alignment/{sample}.srt.bam",
        "genome_alignment/{sample}.srt.filt.bam",
    params:
        sample = lambda wildcards:  wildcards.sample,
        bowtie_path = config["bowtie_path"],
        bowtie_genome_index = config["bowtie_genome_index"],
        samtools_path = config["samtools_path"],
    
    resources: cpus="10", maxtime="2:00:00", mem_mb="60gb",

    shell: """
        {params.bowtie_path} \
            -x {params.bowtie_genome_index} \
            -U {input} -p 12 \
            --very-sensitive-local \
            --un genome_alignment/{params.sample}.unalign.fastq \
            -S genome_alignment/{params.sample}.aln.sam 2>genome_alignment/{params.sample}_genome.log.txt
        # sort and index 
        {params.samtools_path} view -Sb genome_alignment/{params.sample}.aln.sam | \
            {params.samtools_path} sort -@ 4 - > genome_alignment/{params.sample}.srt.bam
        {params.samtools_path} index genome_alignment/{params.sample}.srt.bam
        
        # filter by gap presence 
        {params.samtools_path} view -h genome_alignment/{params.sample}.srt.bam | \
            awk 'BEGIN {{OFS="\t"}} $1 ~ /^@/ || ($0 !~ /XG:i:[^0]/ && $0 !~ /XO:i:[^0]/) {{print $0}}' | \
            {params.samtools_path} view -Sb -o genome_alignment/{params.sample}.srt.filt.bam
        {params.samtools_path} index genome_alignment/{params.sample}.srt.filt.bam
"""



# define function selecting input FASTQ file for alignment
def get_genome_counts_input(wildcards):
    if USE_UMITOOLS:
        return f"genome_alignment/{wildcards.sample}.srt.filt.dedup.bam"
    return f"genome_alignment/{wildcards.sample}.srt.filt.bam"


rule genome_counts:
    input:  
        (expand("genome_alignment/{sample}.srt.filt.dedup.bam", sample=sample_list) 
        if USE_UMITOOLS 
        else expand("genome_alignment/{sample}.srt.filt.bam", sample=sample_list)),
    output: 
        "genome_counts/featurecounts.readcounts.ann.tsv",
        "genome_counts/featurecounts.readcounts_tpm.tsv",
        "genome_counts/featurecounts.readcounts_tpm.ann.tsv",
        "genome_counts/featurecounts.readcounts.raw.tsv.summary"

    params:
        featurecounts = config['featurecounts_path'],
        layout = config["layout"],
        pair_flag = "-p" if config["layout"]=="paired" else "",
        strand = config['featurecounts_strand'],
        gtf = config['annotation_gtf'],
        fc_ann_script = config['featurecounts_annscript'],
    conda:
        "env_config/featurecounts.yaml",

    resources: cpus="10", maxtime="8:00:00", mem_mb="100gb",

    shell: """
        {params.featurecounts} -T 32 -Q 10 {params.pair_flag} -s {params.strand}  -a {params.gtf} -o genome_counts/featurecounts.readcounts.raw.tsv {input}
        sed s/"genome_alignment\/"//g genome_counts/featurecounts.readcounts.raw.tsv| sed s/".srt.bam"//g| tail -n +2 > genome_counts/featurecounts.readcounts.tsv
        python scripts/readcnt_to_rpkmtpm.py genome_counts/featurecounts.readcounts.tsv {params.layout}
        python {params.fc_ann_script} {params.gtf} genome_counts/featurecounts.readcounts.tsv > genome_counts/featurecounts.readcounts.ann.tsv
        python {params.fc_ann_script} {params.gtf} genome_counts/featurecounts.readcounts_tpm.tsv > genome_counts/featurecounts.readcounts_tpm.ann.tsv
"""

rule alignment_metrics_counts:
    input:  
        expand("genome_alignment/{sample}.srt.filt.bam", sample=sample_list),
        expand("genome_alignment/{sample}.srt.filt.dedup.bam", sample=sample_list) if USE_UMITOOLS else [],
        expand("mirbase_alignment/{sample}.srt.bam", sample=sample_list),
        expand("mirbase_alignment/{sample}.srt.dedup.bam", sample=sample_list) if USE_UMITOOLS else [],
        "genome_counts/featurecounts.readcounts.raw.tsv.summary"

    output: 
        "metrics/mirna_genome_alignment_metrics.tsv",
        "metrics/mirna_genome_alignment_metrics.xlsx",

    params:
        use_umi = USE_UMITOOLS,
    conda:
        "env_config/featurecounts.yaml",
    resources: cpus="1", maxtime="8:00:00", mem_mb="2gb",
    shell: """
        mkdir -p metrics
        if [ "{params.use_umi}" = "true" ]; then
            python scripts/qc_metrics-umi.py umi_reads mirbase_alignment genome_alignment
        else
            python scripts/qc_metrics-non-umi.py mirbase_alignment genome_alignment
        fi
"""


rule pca_plots:
    input: "mirbase_counts/mirbase.readcounts.tsv",

    output:
        "plots/PCA_1_vs_2.png",
        "plots/PCA_Variance_Bar_Plot.png",
        "plots/Gene_Variance_Plot.png",

    params:
        num_genes = 500,
        pca_plot_script = config['pca_plot_script'],
        
    conda:
        # uses a subset of the packages that featurecounts does
        "env_config/pcaplot.yaml",

    resources: cpus="1", maxtime="1:00:00", mem_mb="2gb",

    shell: """
        python {params.pca_plot_script} \
        mirbase_counts/mirbase.readcounts.tsv \
        plots \
        --genes_considered {params.num_genes} 
#        --color_file sample_ref/sample_colors_hex.tsv
    """


