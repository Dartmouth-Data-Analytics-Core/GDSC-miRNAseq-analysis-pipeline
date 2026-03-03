#----- Import libraries
import csv

#----- Generate run file which is required for clover-seq rules
def generate_runfile(sample_file):
    with open(sample_file, 'r') as infile, open("runfile.txt", 'w') as outfile:
        reader = csv.DictReader(infile)
        for row in reader:
            sample_id = row["sample_id"]
            group = "dummy"   # <- hardcoded value
            outfile.write(f"{sample_id} {group} contamination\n")

#----- Run function
generate_runfile(config["sample_csv"])

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
    threads: 8
    resources: 
        maxtime="2:00:00", 
        mem_mb="60gb",
    message: "Aligning {wildcards.sample} reads to padded mature miRNA sequences with Bowtie2."
    log: "alignment_logs/mirbase_mature_padded/{sample}.bowtie2.mature.log"
    shell: """

        #----- Make logs subdirectory
        mkdir -p alignment_logs/mirbase_mature_padded

        #----- Run Bowtie1 with unpadded reference
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
            --un mirbase_alignment/{params.sample}.mature.unalign.fastq \
            -S mirbase_alignment/{params.sample}.aln.sam 2> {log}

        #----- Subset reads for aligned length > 16 & < 28bp & any reads with gaps (XO/XG tags)
        {params.samtools_path} view -h mirbase_alignment/{params.sample}.aln.sam | \
            awk 'BEGIN {{OFS="\t"}} $1 ~ /^@/ || ((length($10) > 16 && length($10) <= 28) && ($0 !~ /XG:i:[^0]/ && $0 !~ /XO:i:[^0]/)) {{print $0}}' | \
            samtools view -Sb - > mirbase_alignment/{params.sample}.bam

        #----- Filter for any reads with MAPQ <=1
        {params.samtools_path} view -h -q 2 mirbase_alignment/{params.sample}.bam > mirbase_alignment/{params.sample}.sub.bam
        
        #----- Filter for any reads with > 2 mismatches 
        {params.samtools_path} view -h mirbase_alignment/{params.sample}.sub.bam | \
            awk 'BEGIN {{OFS="\t"}} /^@/ || ($0 ~ /NM:i:[0-2]($|\t)/)' | \
            samtools view -b > mirbase_alignment/{params.sample}.sub2.bam
        
        #----- Sort and index BAM file 
        {params.samtools_path} sort -@ 4 mirbase_alignment/{params.sample}.sub2.bam > mirbase_alignment/{params.sample}.mature.srt.bam
        {params.samtools_path} index mirbase_alignment/{params.sample}.mature.srt.bam

        #----- Remove intermediate bam files 
        rm -rf mirbase_alignment/{params.sample}.bam
        rm -rf mirbase_alignment/{params.sample}.sub.bam
        rm -rf mirbase_alignment/{params.sample}.sub2.bam 

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

#----- Rule to count miRNAs (mature only)
rule mirbase_count:
    """
    Count miRNAs
    """
    input:
        expand("mirbase_alignment/{sample}.mature.srt.bam.idxstats", sample=sample_list),
    output:
        rawCounts = "miRNA_Quant/mature/mature_mirbase.readcounts.tsv",
        tpmCount = "miRNA_Quant/mature/mature_mirbase.readcounts_tpm.tsv",
    resources: 
        cpus="10", 
        maxtime="2:00:00", 
        mem_mb="60gb",
    message: "Counting mature miRNAs."
    shell: """

    #----- Make subdirectory
    mkdir -p miRNA_Quant/mature

    #----- Collate counts
    echo -ne mirbase_ID"\t"Length"\t" > miRNA_Quant/mature/mature_mirbase.readcounts.tsv
    echo {input} | tr " " "\t"| sed s/"mirbase_alignment\/"//g| sed s/".srt.bam.idxstats"//g >> miRNA_Quant/mature/mature_mirbase.readcounts.tsv
    paste {input}| awk -f scripts/mirbase_counts.awk >> miRNA_Quant/mature/mature_mirbase.readcounts.tsv

    #----- Run TPM normalization 
    python scripts/mirbase-readcnt_to_tpm.py miRNA_Quant/mature/mature_mirbase.readcounts.tsv
"""


