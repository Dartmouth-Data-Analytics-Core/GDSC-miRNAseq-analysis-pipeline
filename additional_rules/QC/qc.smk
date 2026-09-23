#----- Rule to make fastqc config
rule make_fastqc_config:
    input:
        samplesheet = config["sample_csv"]
    output:
        "fastQC/fastqc_multiqc_config.yaml"
    params:
        layout = config["layout"]
    container: "docker://ghcr.io/dartmouth-data-analytics-core/fastqc:2.0"
    resources:
        cpus="10", 
        maxtime="2:00:00", 
        mem_mb=61440,
    message: "Building fastQC config."
    shell: """
    
        #----- Build fastqc multiqc config
        echo "sample_names_replace:" > {output}

        #----- Skip header and process lines
        tail -n +2 {input} | while IFS=',' read -r sample fastq1 fastq2; do
            base1=$(basename "$fastq1" .fastq.gz)

            if [[ "{params.layout}" == "paired" ]]; then
                base2=$(basename "$fastq2" .fastq.gz)
                echo "  \\"$base1\\": \\"${{sample}}_Forward\\"" >> {output}
                echo "  \\"$base2\\": \\"${{sample}}_Reverse\\"" >> {output}
            elif [[ "{params.layout}" == "single" ]]; then
                echo "  \\"$base1\\": \\"${{sample}}\\"" >> {output}
            fi
        done
    
    """


#----- Rule to run fastqc on trimmed reads
rule trimmed_fastqc:
    input:
        "trimming/{sample}.R1.trim.fastq.gz"
    output:
        report = "fastQC/{sample}.R1.trim_fastqc.html",
        zipFile = "fastQC/{sample}.R1.trim_fastqc.zip"
    params:
        sample = lambda wildcards:  wildcards.sample,
        fastqc_path = config["fastqc_path"]
    container: "docker://ghcr.io/dartmouth-data-analytics-core/fastqc:2.0"
    threads: 8
    resources:
        cpus="10", 
        maxtime="2:00:00", 
        mem_mb=61440,
    message: "Running {wildcards.sample} trimmed fastQC."
    shell: """
    
        #----- Make output directory
        mkdir -p fastQC

        #----- Run fastQC
        fastqc \
            --threads {threads} \
            -o fastQC \
            {input}
    """
