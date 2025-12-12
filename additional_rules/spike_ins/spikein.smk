rule spikein_bbduk:
    input:
        "trimming/{sample}.R1.trim.fastq.gz"
    output:
        stats = "spikein_alignment/{sample}.stats",
        unmapped = "spikein_alignment/{sample}.unmapped.fastq.gz"
    params:
        spikein_ref = config["spikein_reference_core"]  
    resources: cpus="10", maxtime="4:00:00", mem_mb="60gb",
    threads: 8
    conda:
        "../../env_config/bbmap.yaml"
    shell: """
        mkdir -p spikein_alignment
        
        bbduk.sh \
            in={input} \
            out={output.unmapped} \
            ref={params.spikein_ref} \
            stats={output.stats} \
            k=13 \
            maskmiddle=f \
            rcomp=f \
            hdist=0 \
            edist=0 \
            threads={threads}
    """

rule spikein_counts:
    input:
        expand("spikein_alignment/{sample}.stats", sample=sample_list)
    output:
        "spikein_counts/spikein.readcounts.tsv"
    params:
        samples = " ".join(sample_list)
    resources: cpus="10", maxtime="4:00:00", mem_mb="60gb",
    shell: """
        mkdir -p spikein_counts
        python scripts/parse_bbduk_spikein_stats.py \
            --stats {input} \
            --samples "{params.samples}" \
            --output {output[0]}
        
    """


rule mappingBowtieSpikeIns:
    input:  
        "trimming/{sample}.R1.trim.fastq.gz"
    output: 
        map = "spikein_alignment/{sample}.map",
        unmapped = "spikein_alignment/{sample}.unmapped.bowtie.fastq.gz"
    threads: 12
    resources: cpus="10", maxtime="4:00:00", mem_mb="60gb",
    log:    "spikein_alignment/{sample}.log" 
    conda:
        "../../env_config/bowtie1.yaml"
    shell:
        """
        mkdir -p spikein_alignment
        gunzip -c {input} > spikein_alignment/{wildcards.sample}.fastq
        bowtie --threads {threads} -q -k1 --fullref --best -v0 --norc \
        --un >(gzip > '{output.unmapped}') refs/spikeins/spikeins_full \
        spikein_alignment/{wildcards.sample}.fastq \
        > '{output.map}' 2> {log}
        rm spikein_alignment/{wildcards.sample}.fastq
        """


rule normalize_data_spikein:
    input:
        spikein = "spikein_counts/spikein.readcounts.tsv",
        mirbase_counts = "mirbase_counts/mirbase.readcounts.tsv"
    output:
        "spikein_metrics/spikein_detection_metrics.tsv",
        "spikein_metrics/normalized_scalefactor_mirbase_counts.tsv"
    conda:
        "../../env_config/r_env.yaml"
    resources: cpus="10", maxtime="4:00:00", mem_mb="60gb",
    shell:
        """
        mkdir -p spikein_metrics
        Rscript scripts/normalize_spikein.R
        """
