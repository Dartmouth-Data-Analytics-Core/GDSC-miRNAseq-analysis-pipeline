#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# GDSC miRNA Pipeline v3
#
# Pipeline for the quantification of miRNAs, isomiRs, and other small RNAs
#
# TO DO
# - Add script to collapse isomirs down to their family to provide a decent proxy for mature counts
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
import pandas as pd

#----- set config file
configfile: "config.yaml"
USE_SPIKEINS = config.get("use_spikeins", False)
USE_UMITOOLS = config.get("use_umitools", False)

#----- read in sample data
samples_df = pd.read_csv(config["sample_csv"]).set_index("sample_id", drop=False)
sample_list = list(samples_df['sample_id'])

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# MAIN PIPELINE RULES
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#----- Additional rules based on config file
if USE_SPIKEINS:
    include: "additional_rules/spike_ins/spikein.smk"
if USE_UMITOOLS:
    include: "additional_rules/umitools/umi_extract.smk"

#----- Main pipieline execution
rule all:
    input:
        #== Trimming and UMI Deduplication outputs
        expand("trimming/{sample}.R1.trim.fastq.gz", sample=sample_list),
        expand("trimming/{sample}.cutadapt.report", sample=sample_list),    
        expand("umi_reads/{sample}.umi.fastq.gz", sample=sample_list) if USE_UMITOOLS else [],
        expand("collapsed/{sample}.seqcluster.fastq.gz", sample=sample_list),
        expand("collapsed/{sample}.seqcluster.hairpin.aln.srt.bam", sample=sample_list),

        #== mirbase mature alignment and metrics and hairpin alignment and metrics outputs
        expand("mirbase_alignment/mature/{sample}.mature.srt.bam", sample=sample_list),
        expand("mirbase_alignment/mature/{sample}.mature.srt.dedup.bam", sample=sample_list) if USE_UMITOOLS else [],
        expand("mirbase_alignment/mature/{sample}.mature.unalign.fastq", sample=sample_list),
        expand("mirbase_alignment/mature/{sample}.mature.srt.bam.idxstats", sample=sample_list),
        expand("mirbase_alignment/mature/{sample}.mature.srt.bam.flagstat", sample=sample_list),
        "mirbase_counts/mature_mirbase.readcounts.tsv",
        "mirbase_counts/mature_mirbase.readcounts_tpm.tsv",
        
        #== mirbase hairpin alignments and metrics
        expand("mirbase_alignment/hairpin/{sample}.hairpin.srt.bam", sample=sample_list),
        expand("mirbase_alignment/hairpin/{sample}.hairpin.srt.dedup.bam", sample=sample_list) if USE_UMITOOLS else [],
        expand("mirbase_alignment/hairpin/{sample}.hairpin.unalign.fastq", sample=sample_list),
        expand("mirbase_alignment/hairpin/{sample}.hairpin.srt.bam.idxstats", sample=sample_list),
        expand("mirbase_alignment/hairpin/{sample}.hairpin.srt.bam.flagstat", sample=sample_list),

        #== mirtop
        expand("mirtop/{sample}.hairpin.gff", sample=sample_list),

        #== Genome alignment and metrics outputs
        expand("genome_alignment/{sample}.srt.bam", sample=sample_list),
        expand("genome_alignment/{sample}.srt.filt.bam", sample=sample_list),
        expand("genome_alignment/{sample}.srt.filt.dedup.bam", sample=sample_list) if USE_UMITOOLS else [],
        "genome_counts/featurecounts.readcounts.ann.tsv",
        "genome_counts/featurecounts.readcounts_tpm.tsv",
        "genome_counts/featurecounts.readcounts_tpm.ann.tsv",
        "metrics/mirna_genome_alignment_metrics.tsv",
        
        #== Spike-in, QC, and plot outputs
        expand("spikein_alignment/{sample}.unmapped.bowtie.fastq.gz", sample=sample_list) if USE_SPIKEINS else [],
        "spikein_counts/spikein.readcounts.tsv" if USE_SPIKEINS else [],
        "spikein_metrics/normalized_scalefactor_mirbase_counts.tsv" if USE_SPIKEINS else [],
        "plots/PCA_1_vs_2.png"
    output:
        "multiqc_report.html"          
    conda:
        "env_config/multiqc.yaml",
    resources: cpus="10", maxtime="2:00:00", mem_mb="60gb",
    params:
        multiqc=config["multiqc_path"],
        use_umi = USE_UMITOOLS,
    shell: """
        if [ "{params.use_umi}" = "true" ]; then
            {params.multiqc} -v -c multiqc_config.yaml genome_alignment mirbase_alignment genome_counts mirbase_counts umi_reads
        else
            {params.multiqc} -v -c multiqc_config.yaml genome_alignment mirbase_alignment genome_counts mirbase_counts
        fi
"""

#----- Rule to execute trimming
rule trimming:
    output: 
        "trimming/{sample}.R1.trim.fastq.gz",
        "trimming/{sample}.cutadapt.report"
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
    shell: """
        cutadapt \
            -o trimming/{params.sample}.R1.trim.fastq.gz \
            {params.fastq_file_1} \
            -m 1 \
            {params.nextseq_trim} \
            -j {resources.cpus} \
            -q 30 \
            --max-n 0.8 \
            -a {params.adapter_3prime} \
            --trim-n > trimming/{params.sample}.cutadapt.report
    """

#----- Function selecting input FASTQ file for alignment
def get_alignment_input(wildcards):
    if USE_UMITOOLS:
        return f"umi_reads/{wildcards.sample}.umi.fastq.gz"
    return f"trimming/{wildcards.sample}.R1.trim.fastq.gz"

#----- Rule to collapse reads with seqcluster
rule seqcluster:
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
    shell: """
    
        #----- Run seqcluster to collapse reads
        seqcluster \
            collapse \
            -m 1 \
            --min_size 15 \
            -f {input} \
            -o collapsed
        
        mv collapsed/{params.sample}.R1.trim_trimmed.fastq collapsed/{params.sample}.seqcluster.fastq
        gzip collapsed/{params.sample}.seqcluster.fastq
    
    """

#----- Rule to align collapsed reads to hairpins
rule collapsed_hairpin_aln:
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
    shell: """
    
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
            2>| >(tee collapsed/{params.sample}_aln.out >&2) \
            | {params.samtools_path} view -@ 24 -bS - \
            | {params.samtools_path} sort -@ 24 -o {output.collapsed_aln}

        {params.samtools_path} index {output.collapsed_aln}
    
    """

#----- Rule to align reads to non-padded mature miRNA reference
rule mirbase_mature_aln:
    input:
        get_alignment_input
    output:
        aligned = "mirbase_alignment/mature/{sample}.mature.srt.bam",
        unaligned = "mirbase_alignment/mature/{sample}.mature.unalign.fastq",
    params:
        sample = lambda wildcards:  wildcards.sample,
        bowtie1_path = config["bowtie1_path"],
        bowtie1_index = config["bowtie1_mature_index"],
        samtools_path = config["samtools_path"]
    threads: 8
    resources: 
        maxtime="2:00:00", 
        mem_mb="60gb",
    shell: """

        #----- Run Bowtie1 with unpadded reference
        {params.bowtie1_path} \
            -t \
            -k 50 \
            --best \
            --strata \
            -e 99999 \
            --chunkmbs 2048 \
            --sam \
            -p {threads} \
            -x {params.bowtie1_index} \
            --un {output.unaligned} \
            {input} \
            mirbase_alignment/mature/{params.sample}.mature.sam 2> mirbase_alignment/mature/{params.sample}.mature.bowtie1.log

        #----- Filter for alignment lengths between 16 and 28 bp, convert to bam
        {params.samtools_path} \
            view -h mirbase_alignment/mature/{params.sample}.mature.sam | \
            awk 'BEGIN {{OFS="\t"}} $1 ~ /^@/ || ((length($10) > 16 && length($10) <= 28))' | \
            {params.samtools_path} view -Sb - > mirbase_alignment/mature/{params.sample}.mature.bam

        #----- Sort and index filtered bam file
        {params.samtools_path} \
            sort -@ 4 mirbase_alignment/mature/{params.sample}.mature.bam > {output.aligned}
        
        {params.samtools_path} \
            index {output.aligned}

        #----- Clean temp files
        rm -rf mirbase_alignment/mature/{params.sample}.mature.sam
        rm -rf mirbase_alignment/mature/{params.sample}.mature.bam


"""

#----- Rule to align unaligned reads to hairpin database
rule mirbase_hairpin_aln:
    input:
        "mirbase_alignment/mature/{sample}.mature.unalign.fastq",
    output:
        aligned_hp = "mirbase_alignment/hairpin/{sample}.hairpin.srt.bam",
        unaligned_hp = "mirbase_alignment/hairpin/{sample}.hairpin.unalign.fastq",
    params:
        sample = lambda wildcards:  wildcards.sample,
        bowtie1_path = config["bowtie1_path"],
        bowtie1_hairpin_index = config["bowtie1_hairpin_index"],
        samtools_path = config["samtools_path"],
    threads: 8
    resources: 
        maxtime="2:00:00", 
        mem_mb="60gb",
    shell: """

        #----- Run Bowtie1 with unpadded reference
        {params.bowtie1_path} \
            -t \
            -k 50 \
            --best \
            --strata \
            -e 99999 \
            --chunkmbs 2048 \
            --sam \
            -p {threads} \
            -x {params.bowtie1_hairpin_index} \
            --un {output.unaligned_hp} \
            {input} \
            mirbase_alignment/hairpin/{params.sample}.hairpin.sam 2> mirbase_alignment/hairpin/{params.sample}.hairpin.bowtie1.log

        #----- Filter for alignment lengths between 16 and 28 bp, convert to bam
        {params.samtools_path} \
            view -h mirbase_alignment/hairpin/{params.sample}.hairpin.sam | \
            awk 'BEGIN {{OFS="\t"}} $1 ~ /^@/ || ((length($10) > 16 && length($10) <= 28))' | \
            {params.samtools_path} view -Sb - > mirbase_alignment/hairpin/{params.sample}.hairpin.bam

        #----- Sort and index filtered bam file
        {params.samtools_path} \
            sort -@ 4 mirbase_alignment/hairpin/{params.sample}.hairpin.bam > {output.aligned_hp}
        
        {params.samtools_path} \
            index {output.aligned_hp}

        #----- Clean temp files
        rm -rf mirbase_alignment/hairpin/{params.sample}.hairpin.sam
        rm -rf mirbase_alignment/hairpin/{params.sample}.hairpin.bam

"""

#----- Define function selecting input BAM file for mirbase_stats (mature and hairpin)
def get_mirbase_mature_stats_input(wildcards):
    if USE_UMITOOLS:
        return f"mirbase_alignment/mature/{wildcards.sample}.mature.srt.dedup.bam"
    return f"mirbase_alignment/mature/{wildcards.sample}.mature.srt.bam"
def get_mirbase_hairpin_stats_input(wildcards):
    if USE_UMITOOLS:
        return f"mirbase_alignment/hairpin/{wildcards.sample}.hairpin.srt.dedup.bam"
    return f"mirbase_alignment/hairpin/{wildcards.sample}.hairpin.srt.bam"

#----- Calculate stats for alignments
rule mature_mirbase_stats:
    input: 
        matureStats = get_mirbase_mature_stats_input,
        hairpinStats = get_mirbase_hairpin_stats_input
    output:
        mature_idx = "mirbase_alignment/mature/{sample}.mature.srt.bam.idxstats",
        mature_flagstat = "mirbase_alignment/mature/{sample}.mature.srt.bam.flagstat",
        hairpin_idx = "mirbase_alignment/hairpin/{sample}.hairpin.srt.bam.idxstats",
        hairpin_flagstat = "mirbase_alignment/hairpin/{sample}.hairpin.srt.bam.flagstat",
    params:
        sample = lambda wildcards:  wildcards.sample,
        samtools_path = config["samtools_path"],
    resources: 
        cpus="10", 
        maxtime="2:00:00", 
        mem_mb="60gb",
    shell: """
    {params.samtools_path} idxstats {input.matureStats} > {output.mature_idx}
    {params.samtools_path} flagstat {input.matureStats} > {output.mature_flagstat}
    {params.samtools_path} idxstats {input.hairpinStats} > {output.hairpin_idx}
    {params.samtools_path} flagstat {input.hairpinStats} > {output.hairpin_flagstat}
"""



#----- Rule to run miRtop
rule miRtop:
    input:
        collapsed_aln = "collapsed/{sample}.seqcluster.hairpin.aln.srt.bam"
    output:
        mirtop_gff = "mirtop/{sample}.hairpin.gff"
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
    shell: """

        #----- Run miRtop GFF
        mirtop gff \
            --add-extra \
            --sps {params.sps} \
            --hairpin {params.hairpin_fa} \
            --gtf {params.hairpin_gff} \
            -o mirtop \
            {input.collapsed_aln} &&
        mv mirtop/{params.sample}.seqcluster.hairpin.aln.srt.gff mirtop/{params.sample}.hairpin.gff
        
        #----- Run miRtop counts
        mirtop counts \
            -o mirtop \
            --hairpin {params.hairpin_fa} \
            --gff {output.mirtop_gff} \
            --gtf {params.hairpin_gff}

        #----- Run miRtop stats
        #mirtop stats \
        #    {params.sample}.hairpin.gff \
        #    -o mirtop/stats \
        #    
"""


#----- Rule to count miRNAs (mature only)
rule mirbase_count:
    input:
        expand("mirbase_alignment/mature/{sample}.mature.srt.bam.idxstats", sample=sample_list),
    output:
        rawCounts = "mirbase_counts/mature_mirbase.readcounts.tsv",
        tpmCount = "mirbase_counts/mature_mirbase.readcounts_tpm.tsv",
    params:
    resources: 
        cpus="10", 
        maxtime="2:00:00", 
        mem_mb="60gb",
    shell: """
    echo -ne mirbase_ID"\t"Length"\t" > mirbase_counts/mature_mirbase.readcounts.tsv
    echo {input} | tr " " "\t"| sed s/"mirbase_alignment\/"//g| sed s/".srt.bam.idxstats"//g >> mirbase_counts/mature_mirbase.readcounts.tsv
    paste {input}| awk -f scripts/mirbase_counts.awk >> mirbase_counts/mature_mirbase.readcounts.tsv

    # run TPM normalization 
    python scripts/mirbase-readcnt_to_tpm.py mirbase_counts/mature_mirbase.readcounts.tsv
"""

#----- Rule to conduct alignment to genome
rule genome_alignment:
    input: 
        "mirbase_alignment/hairpin/{sample}.hairpin.unalign.fastq",
    output:
        "genome_alignment/{sample}.srt.bam",
        "genome_alignment/{sample}.srt.filt.bam",
    params:
        sample = lambda wildcards:  wildcards.sample,
        bowtie2_path = config["bowtie2_path"],
        bowtie2_genome_index = config["bowtie2_genome_index"],
        samtools_path = config["samtools_path"],
    
    resources: cpus="10", maxtime="2:00:00", mem_mb="60gb",

    shell: """

        #----- Genome alignment
        {params.bowtie2_path} \
            -x {params.bowtie2_genome_index} \
            -U {input} -p 12 \
            --very-sensitive-local \
            --un genome_alignment/{params.sample}.unalign.fastq \
            -S genome_alignment/{params.sample}.aln.sam 2>genome_alignment/{params.sample}_genome.log.txt
        
        #----- Sort and index 
        {params.samtools_path} view -Sb genome_alignment/{params.sample}.aln.sam | \
            {params.samtools_path} sort -@ 4 - > genome_alignment/{params.sample}.srt.bam
        {params.samtools_path} index genome_alignment/{params.sample}.srt.bam
        
        #----- Filter by gap presence 
        {params.samtools_path} view -h genome_alignment/{params.sample}.srt.bam | \
            awk 'BEGIN {{OFS="\t"}} $1 ~ /^@/ || ($0 !~ /XG:i:[^0]/ && $0 !~ /XO:i:[^0]/) {{print $0}}' | \
            {params.samtools_path} view -Sb -o genome_alignment/{params.sample}.srt.filt.bam
        {params.samtools_path} index genome_alignment/{params.sample}.srt.filt.bam
"""



#----- Define function selecting input FASTQ file for alignment
def get_genome_counts_input(wildcards):
    if USE_UMITOOLS:
        return expand("genome_alignment/{sample}.srt.filt.dedup.bam", sample=sample_list)
    return expand("genome_alignment/{sample}.srt.filt.bam", sample=sample_list)

#----- Count genome alignments
rule genome_counts:
    input:  
        get_genome_counts_input
    output: 
        rawAnno = "genome_counts/featurecounts.readcounts.ann.tsv",
        tpm = "genome_counts/featurecounts.readcounts_tpm.tsv",
        tpmAnno = "genome_counts/featurecounts.readcounts_tpm.ann.tsv",
        summary = "genome_counts/featurecounts.readcounts.raw.tsv.summary"

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

        #----- Count
        {params.featurecounts} \
            -T 32 \
            -Q 10 {params.pair_flag} \
            -s {params.strand} \
            -a {params.gtf} \
            -o genome_counts/featurecounts.readcounts.raw.tsv \
            {input}
        
        #----- Clean
        sed s/"genome_alignment\/"//g genome_counts/featurecounts.readcounts.raw.tsv | \
            sed s/".srt.bam"//g| tail -n +2 > genome_counts/featurecounts.readcounts.tsv
        
        #----- Normalize and annotate
        python scripts/readcnt_to_rpkmtpm.py genome_counts/featurecounts.readcounts.tsv {params.layout}
        python {params.fc_ann_script} {params.gtf} genome_counts/featurecounts.readcounts.tsv > {output.rawAnno}
        python {params.fc_ann_script} {params.gtf} {output.tpm} > {output.tpmAnno}
"""


#----- Rule to get genome alignment metrics
#! TO-DO: need to fix this to work with Bowtie1 logs, for now we turn it off
rule alignment_metrics_counts:
    input:  
        expand("genome_alignment/{sample}.srt.filt.bam", sample=sample_list),
        expand("genome_alignment/{sample}.srt.filt.dedup.bam", sample=sample_list) if USE_UMITOOLS else [],
        expand("mirbase_alignment/mature/{sample}.mature.srt.bam", sample=sample_list),
        expand("mirbase_alignment/mature/{sample}.mature.srt.dedup.bam", sample=sample_list) if USE_UMITOOLS else [],
        "genome_counts/featurecounts.readcounts.raw.tsv.summary"
    output: 
        "metrics/mirna_genome_alignment_metrics.tsv"
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


#----- Rule to plot PCA
rule pca_plots:
    input: "mirbase_counts/mature_mirbase.readcounts.tsv",
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
        mirbase_counts/mature_mirbase.readcounts.tsv \
        plots \
        --genes_considered {params.num_genes} 
#        --color_file sample_ref/sample_colors_hex.tsv
    """


