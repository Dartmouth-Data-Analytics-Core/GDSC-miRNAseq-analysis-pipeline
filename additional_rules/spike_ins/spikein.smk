rule spikein_bbduk:
    """
    Generate spike in data
    """
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
    message: "Generating {wildcards.sample} spike-in data with bbmap."
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
    """
    Get spike-in counts
    """
    input:
        expand("spikein_alignment/{sample}.stats", sample=sample_list)
    output:
        "spikein_counts/spikein.readcounts.tsv"
    params:
        samples=lambda wildcards, input: [path.split("/")[-1].replace(".stats","") for path in input]
    resources: cpus="10", maxtime="4:00:00", mem_mb="60gb",
    message: "Generating spike-in counts."
    shell: """
        echo STATS: {input}
        echo SAMPLES: {params.samples}
        mkdir -p spikein_counts
        python scripts/parse_bbduk_spikein_stats.py \
            --stats {input} \
            --samples {params.samples} \
            --output {output}
    """


rule mappingBowtieSpikeIns:
    """
    Align spike-in data
    """
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
    message: "Aligning {wildcards.sample} spike-in data with Bowtie"
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
    """
    Normalize spike-in data
    """
    input:
        spikein = "spikein_counts/spikein.readcounts.tsv",
        mirbase_counts = "mirbase_counts/mirbase.readcounts.tsv"
    output:
        "spikein_metrics/spikein_detection_metrics.tsv",
        "spikein_metrics/normalized_scalefactor_mirbase_counts.tsv"
    conda:
        "../../env_config/r_env.yaml"
    resources: cpus="10", maxtime="4:00:00", mem_mb="60gb",
    params:
        final_volume = config["sample_with_spikein_finalvolume"]
    message: "Normalizing spike-in data"
    shell:
        """
        mkdir -p spikein_metrics
        Rscript scripts/normalize_spikein.R {params.final_volume}
        """
