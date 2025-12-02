rule umitools:
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

    shell: """
        {params.umitools_path} extract \
            --extract-method=regex \
            --bc-pattern='.+(?P<discard_1>AACTGTAGGCACCATCAAT){{s<=2}}(?P<umi_1>.{{12}})(?P<discard_2>.+)' \
            -I {input} \
            -S umi_reads/{params.sample}.umi.fastq.gz \
            -L umi_reads/{params.sample}.umi.log.txt
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

rule genome_dedup:
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

    shell: """
    {params.umitools_path} dedup --method=unique -I genome_alignment/{params.sample}.srt.bam -S genome_alignment/{params.sample}.srt.dedup.bam
"""
