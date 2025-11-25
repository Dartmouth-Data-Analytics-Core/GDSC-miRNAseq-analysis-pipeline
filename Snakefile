#####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# setup environment
#####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
import pandas as pd

# set config file
configfile: "config.yaml"
USE_SPIKEINS = config.get("use_spikeins", False)

# read in sample data
samples_df = pd.read_table(config["sample_tsv"]).set_index("sample_id", drop=False)
sample_list = list(samples_df['sample_id'])

#####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# define rules
#####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

if USE_SPIKEINS:
    include: "rules/spikein.smk"

rule all:
    input:
        expand("trimming/{sample}.R1.trim.fastq.gz", sample=sample_list),
        expand("trimming/{sample}.cutadapt.report", sample=sample_list),    
        expand("umi_reads/{sample}.umi.fastq.gz", sample=sample_list),
        expand("mirbase_alignment/{sample}.srt.bam", sample=sample_list),
        expand("mirbase_alignment/{sample}.srt.dedup.bam", sample=sample_list),
        expand("genome_alignment/{sample}.srt.bam", sample=sample_list),
        expand("genome_alignment/{sample}.srt.dedup.bam", sample=sample_list),
        expand("genome_alignment/{sample}.srt.dedup.filt.bam", sample=sample_list),
        "metrics/mirna_genome_alignment_metrics.tsv",
        "genome_counts/featurecounts.readcounts.ann.tsv",
        "genome_counts/featurecounts.readcounts_tpm.tsv",
        "genome_counts/featurecounts.readcounts_tpm.ann.tsv",
        "mirbase_counts/mirbase.readcounts.tsv",
        "mirbase_counts/mirbase.readcounts_tpm.tsv",
        expand("spikein_alignment/{sample}.unmapped.bowtie2.fastq.gz", sample=sample_list) if USE_SPIKEINS else [],
        "spikein_counts/spikein.readcounts.tsv" if USE_SPIKEINS else [],
        expand("mirbase_alignment/{sample}.srt.dedup.bam.idxstats", sample=sample_list),
        expand("mirbase_alignment/{sample}.srt.dedup.bam.flagstat", sample=sample_list),
        "plots/PCA_1_vs_2.png"
                
    conda:
        "env_config/multiqc.yaml",
    resources: cpus="10", maxtime="2:00:00", mem_mb="60gb",

    params:
        multiqc=config["multiqc_path"],

    output:
        "multiqc_report.html"

    shell: """
        {params.multiqc}  -c multiqc_config.yaml genome_alignment  mirbase_alignment  genome_counts mirbase_counts  umi_reads
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


def get_input_file(wildcards):
    if USE_SPIKEINS:
        return f"spikein_alignment/{wildcards.sample}.unmapped.bowtie2.fastq.gz"
    return f"trimming/{wildcards.sample}.R1.trim.fastq.gz"


rule umitools:
    input: 
        get_input_file,
    output: 
        "umi_reads/{sample}.umi.fastq.gz",
        "umi_reads/{sample}.umi.log.txt",
    params:
        sample = lambda wildcards:  wildcards.sample,
        umitools_path = config["umitools_path"],
        fastq_file_1 = lambda wildcards: samples_df.loc[wildcards.sample, "fastq_1"],
    resources: cpus="10", maxtime="2:00:00", mem_mb="60gb",

    shell: """
        {params.umitools_path} extract \
            --extract-method=regex \
            --bc-pattern='.+(?P<discard_1>AACTGTAGGCACCATCAAT){{s<=2}}(?P<umi_1>.{{12}})(?P<discard_2>.+)' \
            -I {input} \
            -S umi_reads/{params.sample}.umi.fastq.gz \
            -L umi_reads/{params.sample}.umi.log.txt
"""    



rule mirbase_alignment:
    input:
        "umi_reads/{sample}.umi.fastq.gz"
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


rule mirbase_dedup:
    input: 
        "mirbase_alignment/{sample}.srt.bam",
    output:
        "mirbase_alignment/{sample}.srt.dedup.bam",
    params:
        sample = lambda wildcards:  wildcards.sample,
        bowtie_path = config["bowtie_path"],
        bowtie_index = config["bowtie_index"],
        umitools_path = config["umitools_path"],
        samtools_path = config["samtools_path"],

    resources: cpus="10", maxtime="2:00:00", mem_mb="60gb",

    shell: """
    {params.umitools_path} dedup --method=unique -I mirbase_alignment/{params.sample}.srt.bam -S mirbase_alignment/{params.sample}.srt.dedup.bam
    {params.samtools_path} index mirbase_alignment/{params.sample}.srt.dedup.bam
        
"""

rule mirbase_stats:
    input: 
        "mirbase_alignment/{sample}.srt.bam",
        "mirbase_alignment/{sample}.srt.dedup.bam",
    output:
        "mirbase_alignment/{sample}.srt.dedup.bam.idxstats",
        "mirbase_alignment/{sample}.srt.dedup.bam.flagstat"
    params:
        sample = lambda wildcards:  wildcards.sample,
        bowtie_path = config["bowtie_path"],
        bowtie_index = config["bowtie_index"],
        samtools_path = config["samtools_path"],

    resources: cpus="10", maxtime="2:00:00", mem_mb="60gb",

    shell: """
    {params.samtools_path} idxstats mirbase_alignment/{params.sample}.srt.dedup.bam > mirbase_alignment/{params.sample}.srt.dedup.bam.idxstats
    {params.samtools_path} flagstat mirbase_alignment/{params.sample}.srt.dedup.bam > mirbase_alignment/{params.sample}.srt.dedup.bam.flagstat
        
"""


rule mirbase_count:
    input:
        expand("mirbase_alignment/{sample}.srt.dedup.bam.idxstats", sample=sample_list),

    output:
        "mirbase_counts/mirbase.readcounts.tsv",
        "mirbase_counts/mirbase.readcounts_tpm.tsv",
    params:

    resources: cpus="10", maxtime="2:00:00", mem_mb="60gb",

    shell: """
    echo -ne mirbase_ID"\t"Length"\t" > mirbase_counts/mirbase.readcounts.tsv
    echo {input} | tr " " "\t"| sed s/"mirbase_alignment\/"//g| sed s/".srt.dedup.bam.idxstats"//g >> mirbase_counts/mirbase.readcounts.tsv
    paste {input}| awk -f scripts/mirbase_counts.awk >> mirbase_counts/mirbase.readcounts.tsv

    # run TPM normalization 
    python scripts/mirbase-readcnt_to_tpm.py mirbase_counts/mirbase.readcounts.tsv
"""


rule genome_alignment:
    input: 
        "mirbase_alignment/{sample}.unalign.fastq",
    output:
        "genome_alignment/{sample}.srt.bam",
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
        
"""


rule genome_dedup:
    input: 
        "genome_alignment/{sample}.srt.bam",
    output:
        "genome_alignment/{sample}.srt.dedup.bam",
        "genome_alignment/{sample}.srt.dedup.filt.bam",

    params:
        sample = lambda wildcards:  wildcards.sample,
        bowtie_path = config["bowtie_path"],
        bowtie_index = config["bowtie_index"],
        umitools_path = config["umitools_path"],
        samtools_path = config["samtools_path"],
    
    resources: cpus="10", maxtime="2:00:00", mem_mb="60gb",

    shell: """
    {params.umitools_path} dedup --method=unique -I genome_alignment/{params.sample}.srt.bam -S genome_alignment/{params.sample}.srt.dedup.bam
    # filter by length and gap presence 
    {params.samtools_path} view -h genome_alignment/{params.sample}.srt.dedup.bam | \
        awk 'BEGIN {{OFS="\t"}} $1 ~ /^@/ || ((length($10) > 16 && length($10) <= 28) && ($0 !~ /XG:i:[^0]/ && $0 !~ /XO:i:[^0]/)) {{print $0}}' | \
        {params.samtools_path} view -Sb -> genome_alignment/{params.sample}.srt.dedup.filt.bam
    {params.samtools_path} index genome_alignment/{params.sample}.srt.dedup.filt.bam
"""


rule genome_counts:
    input:  
        expand("genome_alignment/{sample}.srt.dedup.filt.bam", sample=sample_list),

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
        expand("genome_alignment/{sample}.srt.dedup.bam", sample=sample_list),
        expand("mirbase_alignment/{sample}.srt.bam", sample=sample_list),
        expand("mirbase_alignment/{sample}.srt.dedup.bam", sample=sample_list),
        "genome_counts/featurecounts.readcounts.raw.tsv.summary"

    output: 
        "metrics/mirna_genome_alignment_metrics.tsv",
        "metrics/mirna_genome_alignment_metrics.xlsx",

    params:
        gtf = config['annotation_gtf'],
    conda:
        "env_config/featurecounts.yaml",

    resources: cpus="1", maxtime="8:00:00", mem_mb="2gb",

    shell: """
        mkdir -p metrics
        python scripts/qc_metrics.py umi_reads mirbase_alignment genome_alignment
"""


rule pca_plots:
    input: "mirbase_counts/mirbase.readcounts.tsv",

    output:
        #"plots/Heatmap_scaled_"+str(num_genes_compared)+"_features.png",
        # there potentially could be more, but this plot must exist. Make sure -p flag has number at least 2 if specified
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


