#----- Import libraries
import csv

#----- Function selecting input FASTQ file for alignment
def get_alignment_input(wildcards):
    if USE_UMITOOLS:
        return f"umi_reads/{wildcards.sample}.umi.fastq.gz"
    return f"trimming/{wildcards.sample}.R1.trim.fastq.gz"


#----- Rule to align reads to padded mature miRNA reference
rule mirbase_padded_aln:
    """
    Align sequences to mirbase mature padded reference
    """
    input:
        get_alignment_input
    output:
        aligned = "mirbase_alignment/{sample}.mature.srt.bam",
        unaligned = "mirbase_alignment/{sample}.mature.unalign.fastq",
    params:
        sample = lambda wildcards:  wildcards.sample,
        bowtie2_path = config["bowtie2_path"],
        padded_mature_index = config["padded_mature_index"],
        samtools_path = config["samtools_path"]
    threads: 12
    resources: 
        maxtime="2:00:00", 
        mem_mb="60gb",
    message: "Aligning {wildcards.sample} reads to padded mature miRNA sequences with Bowtie2."
    log: "alignment_logs/mirbase_mature_padded/{sample}.bowtie2.mature.log"
    shell: """

        #----- Make logs subdirectory
        mkdir -p alignment_logs/mirbase_mature_padded

        #----- Align with Bowtie2 to padded reference
        {params.bowtie2_path} \
            -x {params.padded_mature_index} \
            -U {input} \
            -p {threads} \
            --norc \
            -D 20 \
            -R 3 \
            -N 1 \
            -L 12 \
            -i S,1,0.50 \
            --un {output.unaligned} \
            -S mirbase_alignment/{params.sample}.mature.aln.sam 2> {log}

        #----- Subset reads for aligned length > 16 & < 28bp & any reads with gaps (XO/XG tags)
        {params.samtools_path} \
            view -h mirbase_alignment/{params.sample}.mature.aln.sam | \
            awk 'BEGIN {{OFS="\t"}} $1 ~ /^@/ || ((length($10) > 16 && length($10) <= 28) && ($0 !~ /XG:i:[^0]/ && $0 !~ /XO:i:[^0]/)) {{print $0}}' | \
            samtools view -Sb -o mirbase_alignment/{params.sample}.mature.bam
        
        #----- Filter for any reads with MAPQ <=1
        {params.samtools_path} \
            view -h -q 2 mirbase_alignment/{params.sample}.mature.bam > mirbase_alignment/{params.sample}.mature.sub.bam
        
        #----- Filter for any reads with > 2 mismatches 
        {params.samtools_path} \
            view -h mirbase_alignment/{params.sample}.mature.sub.bam | \
            awk 'BEGIN {{OFS="\t"}} /^@/ || ($0 ~ /NM:i:[0-2]($|\t)/)' | \
            samtools view -b > mirbase_alignment/{params.sample}.mature.sub2.bam
        
        #----- Sort and index BAM file 
        {params.samtools_path} sort -@ 4 mirbase_alignment/{params.sample}.mature.sub2.bam > {output.aligned}
        {params.samtools_path} index {output.aligned}

        # remove intermediate bam files 
        rm -rf mirbase_alignment/{params.sample}.mature.bam
        rm -rf mirbase_alignment/{params.sample}.mature.sub.bam
        rm -rf mirbase_alignment/{params.sample}.mature.sub2.bam 

"""

#----- Define function selecting input BAM file for mirbase_stats (mature and hairpin)
def get_mirbase_mature_stats_input(wildcards):
    if USE_UMITOOLS:
        return f"mirbase_alignment/{wildcards.sample}.mature.srt.dedup.bam"
    return f"mirbase_alignment/{wildcards.sample}.mature.srt.bam"

#----- Calculate stats for alignments
rule mature_mirbase_stats:
    """
    Collate stats for mirbase alignments
    """
    input: 
        matureStats = get_mirbase_mature_stats_input,
    output:
        mature_idx = "mirbase_alignment/{sample}.mature.srt.bam.idxstats",
        mature_flagstat = "mirbase_alignment/{sample}.mature.srt.bam.flagstat",
    params:
        sample = lambda wildcards:  wildcards.sample,
        samtools_path = config["samtools_path"],
    resources: 
        cpus="10", 
        maxtime="2:00:00", 
        mem_mb="60gb",
    message: "Collating {wildcards.sample} mirbase stats with Samtools"
    shell: """
    {params.samtools_path} idxstats {input.matureStats} > {output.mature_idx}
    {params.samtools_path} flagstat {input.matureStats} > {output.mature_flagstat}
"""

#----- Rule to map unaligned to genome
rule genome_alignment:
    """
    Map reads that did not align to padded mature to genome
    """
    input:
        unaligned = "mirbase_alignment/{sample}.mature.unalign.fastq"
    output:
        genomeAln = "genome_alignment/{sample}.genome.srt.filt.bam"
    params:
        sample = lambda wildcards:  wildcards.sample,
        bowtie2_path = config["bowtie2_path"],
        bowtie2_genome_index = config["bowtie2_genome_index"],
        samtools_path = config["samtools_path"]
    threads: 12
    resources:
        maxtime="2:00:00", 
        mem_mb="60gb",
    message: "Aligning {wildcards.sample} mature unaligned reads to genome."
    log: "alignment_logs/genome_alignment/{sample}.bowtie2.genome.log"
    shell: """
    
        #----- Align mature unaligned reads to genome with bowtie2
        {params.bowtie2_path} \
            -x {params.bowtie2_genome_index} \
            -U {input.unaligned} \
            -p {threads} \
            --very-sensitive \
            -S genome_alignment/{params.sample}.genome.aln.sam 2> {log}
        
        #----- Convert to bam
        {params.samtools_path} \
            view -Sb genome_alignment/{params.sample}.genome.aln.sam | \
            {params.samtools_path} sort -@ 4 -o genome_alignment/{params.sample}.genome.srt.bam

        #----- Filter
        {params.samtools_path} \
            view -h genome_alignment/{params.sample}.genome.srt.bam | \
            awk 'BEGIN {{OFS="\t"}} $1 ~ /^@/ || ($0 !~ /XG:i:[^0]/ && $0 !~ /XO:i:[^0]/) {{print $0}}' | \
            {params.samtools_path} view -Sb -o {output.genomeAln}

        #----- Index
        {params.samtools_path} index {output.genomeAln}
    """

#----- Define function selecting input BAM file for mirbase_stats (mature and hairpin)
def get_genome_stats_input(wildcards):
    if USE_UMITOOLS:
        return f"genome_alignment/{wildcards.sample}.genome.srt.filt.dedup.bam"
    return f"genome_alignment/{wildcards.sample}.genome.srt.filt.bam"

#----- Rule to get genome stats
rule genome_stats:
    """
    Collate stats for genome alignments
    """
    input: 
        genomeStats = get_genome_stats_input,
    output:
        genome_idx = "genome_alignment/{sample}.genome.srt.filt.bam.idxstats",
        genome_flagstat = "mirbase_alignment/{sample}.genome.srt.filt.bam.flagstat",
    params:
        sample = lambda wildcards:  wildcards.sample,
        samtools_path = config["samtools_path"],
    resources: 
        cpus="10", 
        maxtime="2:00:00", 
        mem_mb="60gb",
    message: "Collating {wildcards.sample} mirbase stats with Samtools"
    shell: """
    {params.samtools_path} idxstats {input.genomeStats} > {output.genome_idx}
    {params.samtools_path} flagstat {input.genomeStats} > {output.genome_flagstat}
"""

    
#----- Define function selecting input BAM file for featurecounts
def get_genome_featureCounts_input(wildcards):
    if USE_UMITOOLS:
        return expand("genome_alignment/{sample}.genome.srt.filt.dedup.bam", sample=sample_list)
    return expand("genome_alignment/{sample}.genome.srt.filt.bam", sample=sample_list)

#----- Rule to run featurecounts on genome alignment
rule genome_featureCounts:
    """
    Assign genome reads to features
    """
    input:
        genAln = get_genome_featureCounts_input
    output:
        rawCounts = "genome_counts/featurecounts.tsv",
        cleanCounts = "genome_counts/featurecounts.readcounts.tsv",
        biotype = "genome_counts/featurecounts.readcounts.biotype.tsv"
    params:
        featurecounts_path = config["featurecounts_path"],
        layout = config["layout"],
        pair_flag = "-p" if config["layout"]=="paired" else "",
        featurecounts_strand = config["featurecounts_strand"],
        annotation_gtf = config["annotation_gtf"]
    threads: 32
    resources:
        cpus = "10",
        maxtime = "8:00:00",
        mem_mb = "100gb"
    message: "Running featurecounts on genome alignments."
    shell: """
    
        #----- Run FeatureCounts
        {params.featurecounts_path} \
            -T {threads} \
            -Q 10 \
            {params.pair_flag} \
            -s 0 \
            -a {params.annotation_gtf} \
            -o {output.rawCounts} \
            {input.genAln}
        
        #----- Combine counts
        sed s/"genome_alignment\/"//g {output.rawCounts}| sed s/".genome.srt.filt.bam"//g| tail -n +2 > {output.cleanCounts}
        
        #----- Add annotation
        python scripts/add_biotype.py \
            {output.cleanCounts} \
            {params.annotation_gtf} \
            {output.biotype}

        #----- Normalize
        python scripts/readcnt_to_rpkmtpm.py \
            {output.biotype} \
            {params.layout}

    
    """

#----- Rule to get alignment metrics
rule alignment_metrics_counts:
    """
    Get alignment metrics
    """
    input:  
        expand("genome_alignment/{sample}.genome.srt.filt.bam", sample=sample_list),
        expand("genome_alignment/{sample}.genome.srt.filt.dedup.bam", sample=sample_list) if USE_UMITOOLS else [],
        expand("mirbase_alignment/{sample}.mature.srt.bam", sample=sample_list),
        expand("mirbase_alignment/{sample}.mature.srt.dedup.bam", sample=sample_list) if USE_UMITOOLS else [],
        biotypes = "genome_counts/featurecounts.readcounts.biotype.tsv"
    output: 
        "metrics/mirna_genome_alignment_metrics.tsv",
        "metrics/mirna_genome_alignment_metrics.xlsx",
    conda:
        "../../env_config/featurecounts.yaml",
    params:
        use_umi = USE_UMITOOLS,
    resources: cpus="1", maxtime="8:00:00", mem_mb="2gb",
    message: "Calculating metrics"
    shell: """

        #----- Create directory
        mkdir -p metrics

        #----- Calculate mapping metrics
        if [ "{params.use_umi}" = "true" ]; then
            python scripts/qc_metrics-umi.py \
                umi_reads \
                alignment_logs/mirbase_mature_padded \
                mirbase_alignment \
                alignment_logs/genome_alignment
        else
            python scripts/qc_metrics-non-umi.py \
                alignment_logs/mirbase_mature_padded \
                mirbase_alignment \
                alignment_logs/genome_alignment \
                alignment_logs/seqcluster 
        fi

        #----- Calculate biotype metrics
        python scripts/gene_biotype_metrics.py \
            {input.biotypes}
"""


