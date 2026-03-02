# define function selecting input FASTQ file for UMItools 
def get_umitools_input(wildcards):
    if USE_SPIKEINS:
        return f"spikein_alignment/{wildcards.sample}.unmapped.bowtie.fastq.gz"
    return f"trimming/{wildcards.sample}.R1.trim.fastq.gz"

rule umitools:
    """
    UMI extraction
    """
    input: 
        get_umitools_input,
    output: 
        "umi_reads/{sample}.umi.fastq.gz",
        "umi_reads/{sample}.umi.log.txt",
    params:
        sample = lambda wildcards:  wildcards.sample,
        umitools_path = config["umitools_path"],
        fastq_file_1 = lambda wildcards: samples_df.loc[wildcards.sample, "fastq_1"],
    resources: cpus="10", maxtime="2:00:00", mem_mb="60gb",
    message: "Extracting {wildcards.sample} UMIs."
    shell: """
        {params.umitools_path} extract \
            --extract-method=regex \
            --bc-pattern='.+(?P<discard_1>AACTGTAGGCACCATCAAT){{s<=2}}(?P<umi_1>.{{12}})(?P<discard_2>.+)' \
            -I {input} \
            -S umi_reads/{params.sample}.umi.fastq.gz \
            -L umi_reads/{params.sample}.umi.log.txt
"""    

#----- Rule to deduplicate
rule mirbase_dedup:
    """
    Deduplicate mirBase reads
    """
    input: 
        mature = "mirbase_alignment/mature/{sample}.mature.srt.bam",
        hairpin = "mirbase_alignment/hairpin/{sample}.hairpin.srt.bam"
    output:
        mature_dedup = "mirbase_alignment/mature/{sample}.mature.srt.dedup.bam",
        hairpin_dedup = "mirbase_alignment/hairpin/{sample}.hairpin.srt.dedup.bam"
    params:
        sample = lambda wildcards:  wildcards.sample,
        bowtie1_path = config["bowtie1_path"],
        bowtie1_index = config["bowtie1_index"],
        umitools_path = config["umitools_path"],
        samtools_path = config["samtools_path"],
    resources: cpus="10", maxtime="2:00:00", mem_mb="60gb",
    message: "Deduplicating {wildcards.sample} mirbase reads."
    shell: """

    #----- Deduplicate mature and index
    {params.umitools_path} \
        dedup \
        --method=unique \
        -I {input.mature} \
        -S {output.hairpin_dedup}

    {params.samtools_path} \
        index mirbase_alignment/mature/{sample}.mature.srt.dedup.bam

    #----- Deduplicate hairpin and index
    {params.umitools_path} \
        dedup \
        --method=unique \
        -I {input.hairpin} \
        -S {output.hairpin_dedup}

    {params.samtools_path} \
        index mirbase_alignment/hairpin/{sample}.hairpin.srt.dedup.bam
        
"""

rule genome_dedup:
    """
    Deduplicate genome
    """
    input: 
        "genome_alignment/{sample}.srt.filt.bam",
    output:
        "genome_alignment/{sample}.srt.filt.dedup.bam",
    params:
        sample = lambda wildcards:  wildcards.sample,
        bowtie_path = config["bowtie_path"],
        bowtie_index = config["bowtie_index"],
        umitools_path = config["umitools_path"],
        samtools_path = config["samtools_path"],
    resources: cpus="10", maxtime="2:00:00", mem_mb="60gb",
    message: "Deduplicating {wildcards.sample} genome reads."
    shell: """
    {params.umitools_path} dedup --method=unique -I genome_alignment/{params.sample}.srt.filt.bam -S genome_alignment/{params.sample}.srt.filt.dedup.bam
"""
