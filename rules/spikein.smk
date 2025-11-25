rule spikein_bbduk:
    input:
        "trimming/{sample}.R1.trim.fastq.gz"
    output:
        stats = "spikein_alignment/{sample}.stats",
        unmapped = "spikein_alignment/{sample}.unmapped.fastq.gz"
    params:
        spikein_ref = config["spikein_reference_core"]  
    threads: 8
    conda:
        "../env_config/bbmap.yaml"
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
    shell: """
        mkdir -p spikein_counts
        python scripts/parse_bbduk_spikein_stats.py \
            --stats {input} \
            --samples "{params.samples}" \
            --output {output[0]}
        
    """


rule spikein_bowtie_full:
    input: 
        "trimming/{sample}.R1.trim.fastq.gz"
    output: 
        unmapped = "spikein_alignment/{sample}.unmapped.bowtie2.fastq.gz",
    params:
        bowtie_path = config["bowtie_path"],
        bowtie_spikein_index = config["bowtie_spikein_index"]
    threads: 12
    shell: """
        mkdir -p spikein_alignment
        {params.bowtie_path} \
            --threads {threads} \
            -U <(gunzip -c {input}) \
            -k 1 --norc -N 0 -L 13 \
            --un >(gzip -c > {output.unmapped}) \
            -x {params.bowtie_spikein_index} \
            > /dev/null \
            2> spikein_alignment/{wildcards.sample}.bowtie2.log
    """



rule normalize_data_spikein:
    input:
        spikein = "spikein_counts/spikein.readcounts.tsv",
        mirbase_counts = "mirbase_counts/mirbase.readcounts.tsv"
    output:
        "spikein_metrics/spikein_detection_metrics.tsv",
        "spikein_metrics/normalized_scalefactor_mirbase_counts.tsv"
    shell:
        """
        mkdir -p spikein_metrics
        Rscript scripts/normalize_spikein.R
        """