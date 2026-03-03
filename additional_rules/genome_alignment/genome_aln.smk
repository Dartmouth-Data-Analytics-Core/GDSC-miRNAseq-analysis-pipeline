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

#----- Rule to integrate isomir counts and mature counts
rule integrate_mature_and_isomirs:
    """
    Intgrate mature and isomiR counts
    """
    input:
        mature = "miRNA_Quant/mature/mature_mirbase.readcounts.tsv",
        isomirs = "miRNA_Quant/isomiRs/isomirs_collapsed_by_miRNA_counts.csv"
    output:
        integrated = "miRNA_Quant/merged/merged_mature_and_isomir_counts.csv"
    conda: "../../env_config/r_env.yaml"
    resources:
        cpus="10", 
        maxtime="2:00:00", 
        mem_mb="60gb",
    message: "Integrating mature and isomir counts."
    shell: """

        #----- Make subdirectory
        mkdir -p miRNA_Quant/merged

        #----- Integrated isomir counts and mature counts
        Rscript scripts/merge_mature_and_isomirs.R \
            miRNA_Quant/isomiRs/isomirs_collapsed_by_miRNA_counts.csv \
            miRNA_Quant/mature/mature_mirbase.readcounts.tsv \
            miRNA_Quant/merged/merged_mature_and_isomir_counts.csv
    
    """

#----- Filtering
rule tRNA_mapping:
    """
    Aligning reads that did not align to mature miRNA sequences to tRNA/smRNA database (Clover-Seq)
    """
    input:
        unaligned = "mirbase_alignment/{sample}.mature.unalign.fastq"
    output:
        srtBam = "contamination/{sample}.contam.srt.bam",
        unaligned_contam = "contamination/{sample}.contam.unaligned.fastq"
    params:
        sample = lambda wildcards:  wildcards.sample,
        bowtie2_path = config["bowtie2_path"],
        tRNA_database = config["tRNA_database"],
        tRNA_bowtie2_index = config["tRNA_bowtie2_index"],
        samtools_path = config["samtools_path"],
    resources: 
        cpus="10", 
        maxtime="6:00:00", 
        mem_mb="60gb"
    message: "Mapping {wildcards.sample} mature miRNA unmapped reads to tRNA database."
    log: "alignment_logs/clover-seq/{sample}.bowtie2.tRNA.log"
    shell: """
    
        #----- Make log subdirectory
        mkdir -p alignment_logs/clover-seq

        #----- Run Bowtie2
        {params.bowtie2_path} \
            -x {params.tRNA_bowtie2_index} \
            -U {input.unaligned} \
            -D 20 \
            -R 3 \
            -N 1 \
            -L 12 \
            -i S,1,0.50 \
            -k 100 \
            --very-sensitive \
            --np 5 \
            --ignore-qual \
            --un {output.unaligned_contam} \
            -p {resources.cpus} \
            -S contamination/{params.sample}.contam.aln.sam 2> {log}

        #----- subset reads for aligned length > 15 & < 90bp 
        {params.samtools_path} view -h contamination/{params.sample}.contam.aln.sam | \
            awk 'BEGIN {{OFS="\t"}} $1 ~ /^@/ || ((length($10) > 15 && length($10) <= 90))' | \
            {params.samtools_path} view -Sb - > contamination/{params.sample}.contam.bam
        
        #----- filter for any reads with MAPQ <=1
        {params.samtools_path} view -h -q 2 contamination/{params.sample}.contam.bam > contamination/{params.sample}.sub.bam

        #----- Sort and filter the bam file
        {params.samtools_path} sort -@ 4 contamination/{params.sample}.sub.bam > {output.srtBam}
        {params.samtools_path} index {output.srtBam}

        #----- Remove temp files
        rm -rf contamination/{params.sample}.contam.aln.sam
        rm -rf contamination/{params.sample}.contam.bam
        rm -rf contamination/{params.sample}.sub.bam
    """

#----- Rule to quantify contamination smRNAs
rule contamination_count:
    """
    Count smRNA reads (contamination) with Clover-Seq
    """
    input:
        expand("contamination/{sample}.contam.srt.bam", sample=sample_list)
    output:
        groupCounts = "contamination/counts/smRNA_raw_counts_by_group.txt",
        counts = "contamination/counts/smRNA_raw_counts_by_sample.txt",
        subGroupFile = "contamination/counts/subroup_counts.txt"
    conda: "clover-seq"
    params:
        smRNA_count = "scripts/clover-seq/count_all_smRNA.py",
        tRNA_database = config["tRNA_database"]
    resources: 
        cpus="12", 
        maxtime="6:00:00", 
        mem_mb="60gb"
    message: "Counting smRNA contamination."
    shell: """
    
        #----- Run the code to count all tRNA + smRNA
        python {params.smRNA_count} \
            --samplefile=runfile.txt \
            --trnatable={params.tRNA_database}/db-trnatable.txt \
            --ensemblgtf={params.tRNA_database}/genes.gtf \
            --trnaloci={params.tRNA_database}/db-trnaloci.bed \
            --maturetrnas={params.tRNA_database}/db-maturetRNAs.bed \
            --realcountfile={output.counts} \
            --countfile={output.groupCounts} \
            --mismatchfile={output.subGroupFile}    
    """

#!!! DO WE WANT TO ADD SOMETHING HERE TO NOT REPORT THE MIRNAS OUTPUT BY CLOVER-SEQ?
#rule map_piRNA
    
#----- Rule to conduct alignment to genome
rule genome_alignment:
    """
    Aligning all unaligned to genome.
    """
    input: 
        "mirbase_alignment/{sample}.mature.unalign.fastq",
    output:
        "genome_alignment/{sample}.srt.bam",
        "genome_alignment/{sample}.srt.filt.bam",
    params:
        sample = lambda wildcards:  wildcards.sample,
        bowtie2_path = config["bowtie2_path"],
        bowtie2_genome_index = config["bowtie2_genome_index"],
        samtools_path = config["samtools_path"],
    resources: cpus="10", maxtime="2:00:00", mem_mb="60gb",
    message: "Aligning {wildcards.sample} padded mature unaligned reads to full genome."
    log: "alignment_logs/genome_alignment/{sample}.bowtie2.genome.log"
    shell: """

        #----- Make logs subdirectory
        mkdir -p alignment_logs/genome_alignment

        #----- Genome alignment
        {params.bowtie2_path} \
            -x {params.bowtie2_genome_index} \
            -U {input} \
            -p 12 \
            --very-sensitive-local \
            --un genome_alignment/{params.sample}.unalign.fastq \
            -S genome_alignment/{params.sample}.aln.sam 2> {log}
        
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
    """
    Get genome counts
    """
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
        "../../env_config/featurecounts.yaml",
    resources: cpus="10", maxtime="8:00:00", mem_mb="100gb",
    message: "Getting genome counts."
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

"""
#----- Rule to get genome alignment metrics
rule alignment_metrics_counts:
    input:  
        expand("genome_alignment/{sample}.srt.filt.bam", sample=sample_list),
        expand("genome_alignment/{sample}.srt.filt.dedup.bam", sample=sample_list) if USE_UMITOOLS else [],
        expand("mirbase_alignment/{sample}.mature.srt.bam", sample=sample_list),
        expand("mirbase_alignment/{sample}.mature.srt.dedup.bam", sample=sample_list) if USE_UMITOOLS else [],
        "genome_counts/featurecounts.readcounts.raw.tsv.summary"
    output: 
        "metrics/mirna_genome_alignment_metrics.tsv"
    params:
        use_umi = USE_UMITOOLS,
    conda:
        "../../env_config/featurecounts.yaml",
    resources: cpus="1", maxtime="8:00:00", mem_mb="2gb",
    message: "Collating alignment metrics."
    shell: 
        mkdir -p metrics

        if [ "{params.use_umi}" = "true" ]; then
            python scripts/qc_metrics-umi.py umi_reads alignment_logs/mirbase_mature_padded/{sample}.bowtie2.mature.log alignment_logs/genome_alignment/{sample}.bowtie2.genome.log
        else
            python scripts/qc_metrics-non-umi.py alignment_logs/mirbase_mature_padded/{sample}.bowtie2.mature.log alignment_logs/genome_alignment/{sample}.bowtie2.genome.log
        fi

"""
#----- Rule to plot PCA
rule pca_plots:
    """
    Principal component analysis
    """
    input: "miRNA_Quant/mature/mature_mirbase.readcounts.tsv",
    output:
        "plots/PCA_1_vs_2.png",
        "plots/PCA_Variance_Bar_Plot.png",
        "plots/Gene_Variance_Plot.png",
    params:
        num_genes = 500,
        pca_plot_script = config['pca_plot_script'],   
    conda:
        "../../env_config/pcaplot.yaml",
    resources: cpus="1", maxtime="1:00:00", mem_mb="2gb",
    message: "Calculating principal components."
    shell: """
        python {params.pca_plot_script} \
        {input} \
        plots \
        --genes_considered {params.num_genes} 
    """