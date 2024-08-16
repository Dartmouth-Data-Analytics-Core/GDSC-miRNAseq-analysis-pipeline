#####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# setup environment
#####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
import pandas as pd

# set config file
configfile: "config.yaml"

# read in sample data
samples_df = pd.read_table(config["sample_tsv"]).set_index("sample_id", drop=False)
sample_list = list(samples_df['sample_id'])

#####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# define rules
#####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

rule all:
    input:
        expand("umi_reads/{sample}.umi.fastq.gz", sample=sample_list),
        expand("mirbase_alignment/{sample}.srt.bam", sample=sample_list),
        expand("mirbase_alignment/{sample}.srt.dedup.bam", sample=sample_list),
        expand("genome_alignment/{sample}.srt.bam", sample=sample_list),
        expand("genome_alignment/{sample}.srt.dedup.bam", sample=sample_list),
        "metrics/mirna_genome_alignment_metrics.tsv",
        "genome_counts/featurecounts.readcounts.ann.tsv",
        "mirbase_counts/mirbase.readcounts.tsv",
        "plots/PCA_Variance_Bar_Plot.png",
        #"featurecounts/featurecounts.readcounts_fpkm.ann.tsv",
        expand("mirbase_alignment/{sample}.srt.dedup.bam.idxstats", sample=sample_list),
        expand("mirbase_alignment/{sample}.srt.dedup.bam.flagstat", sample=sample_list)
                
    conda:
        "env_config/multiqc.yaml",
    resources: cpus="10", maxtime="2:00:00", mem_mb="60gb",

    params:
        multiqc=config["multiqc_path"],

    output:
        "multiqc_report.html"

    shell: """
        {params.multiqc}  genome_alignment  mirbase_alignment  genome_counts mirbase_counts  umi_reads


"""



rule umitools:
    output: 
        "umi_reads/{sample}.umi.fastq.gz",
        "umi_reads/{sample}.umi.log.txt",
    params:
        sample = lambda wildcards:  wildcards.sample,
        umitools_path = config["umitools_path"],
        fastq_file_1 = lambda wildcards: samples_df.loc[wildcards.sample, "fastq_1"],
        layout=config["layout"],
        #umi_bc_pattern=config["umi_bc_pattern"]

    resources: cpus="10", maxtime="2:00:00", mem_mb="60gb",

    shell: """
        {params.umitools_path} extract --extract-method=regex --bc-pattern='.+(?P<discard_1>AACTGTAGGCACCATCAAT){{s<=2}}(?P<umi_1>.{{12}})(?P<discard_2>.+)' -I {params.fastq_file_1} -S umi_reads/{params.sample}.umi.fastq.gz -L umi_reads/{params.sample}.umi.log.txt
"""


rule mirbase_alignment:
    input: 
        "umi_reads/{sample}.umi.fastq.gz",
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
        {params.bowtie_path}  -x {params.bowtie_index} -U {input}  -p 12  --norc --very-sensitive-local --un mirbase_alignment/{params.sample}.unalign.fastq -S mirbase_alignment/{params.sample}.aln.sam 2>mirbase_alignment/{params.sample}_mirbase.log.txt

        {params.samtools_path} view -Sb mirbase_alignment/{params.sample}.aln.sam | {params.samtools_path} sort -@ 4 - > mirbase_alignment/{params.sample}.srt.bam
        {params.samtools_path} index mirbase_alignment/{params.sample}.srt.bam
        
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
    params:

    resources: cpus="10", maxtime="2:00:00", mem_mb="60gb",

    shell: """
    echo -ne mirbase_ID"\t"Length"\t" > mirbase_counts/mirbase.readcounts.tsv
    echo {input} | tr " " "\t"| sed s/"mirbase_alignment\/"//g| sed s/".srt.dedup.bam.idxstats"//g >> mirbase_counts/mirbase.readcounts.tsv
    paste {input}| awk -f scripts/mirbase_counts.awk >> mirbase_counts/mirbase.readcounts.tsv
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

        {params.samtools_path} view -Sb genome_alignment/{params.sample}.aln.sam | {params.samtools_path} sort -@ 4 - > genome_alignment/{params.sample}.srt.bam
        {params.samtools_path} index genome_alignment/{params.sample}.srt.bam
        
"""


rule genome_dedup:
    input: 
        "genome_alignment/{sample}.srt.bam",
    output:
        "genome_alignment/{sample}.srt.dedup.bam",
    params:
        sample = lambda wildcards:  wildcards.sample,
        bowtie_path = config["bowtie_path"],
        bowtie_index = config["bowtie_index"],
        umitools_path = config["umitools_path"],
        samtools_path = config["samtools_path"],
    
    resources: cpus="10", maxtime="2:00:00", mem_mb="60gb",

    shell: """
    {params.umitools_path} dedup --method=unique -I genome_alignment/{params.sample}.srt.bam -S genome_alignment/{params.sample}.srt.dedup.bam
    {params.samtools_path} index genome_alignment/{params.sample}.srt.dedup.bam
        
"""


rule genome_counts:
    input:  
        expand("genome_alignment/{sample}.srt.bam", sample=sample_list),

    output: 
        "genome_counts/featurecounts.readcounts.ann.tsv",


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
        {params.featurecounts} -T 32 {params.pair_flag} -s {params.strand}  -a {params.gtf} -o genome_counts/featurecounts.readcounts.raw.tsv {input}
        sed s/"genome_alignment\/"//g genome_counts/featurecounts.readcounts.raw.tsv| sed s/".srt.bam"//g| tail -n +2 > genome_counts/featurecounts.readcounts.tsv
        python {params.fc_ann_script} {params.gtf} genome_counts/featurecounts.readcounts.tsv > genome_counts/featurecounts.readcounts.ann.tsv
"""

rule alignment_metrics_counts:
    input:  
        expand("genome_alignment/{sample}.srt.bam", sample=sample_list),

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
        python scripts/qc_metrics.py umi_reads mirbase_alignment genome_alignment > metrics/mirna_genome_alignment_metrics.tsv
        python scripts/qc_metrics_xlsx.py metrics/mirna_genome_alignment_metrics.tsv metrics/mirna_genome_alignment_metrics.xlsx
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
        --genes_considered {params.num_genes} \
#        --color_file sample_ref/sample_colors_hex.tsv
    """


