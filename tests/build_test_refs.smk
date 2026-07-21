rule all:
    input:
        "tests/fixtures/refs/hairpin_subset.1.ebwt",
        "tests/fixtures/refs/mature_subset.1.ebwt",
        "tests/fixtures/refs/tiny_genome.1.bt2",

rule bowtie1_hairpin:
    input: "tests/fixtures/refs/hairpin_subset.fa"
    output: "tests/fixtures/refs/hairpin_subset.1.ebwt"
    conda: "../env_config/bowtie1.yaml"
    container: "docker://ghcr.io/dartmouth-data-analytics-core/bowtie1:2.0"
    shell: "bowtie-build {input} tests/fixtures/refs/hairpin_subset"

rule bowtie1_mature:
    input: "tests/fixtures/refs/mature_subset.fa"
    output: "tests/fixtures/refs/mature_subset.1.ebwt"
    conda: "../env_config/bowtie1.yaml"
    container: "docker://ghcr.io/dartmouth-data-analytics-core/bowtie1:2.0"
    shell: "bowtie-build {input} tests/fixtures/refs/mature_subset"

rule bowtie2_genome:
    input: "tests/fixtures/refs/tiny_genome.fa"
    output: "tests/fixtures/refs/tiny_genome.1.bt2"
    conda: "../env_config/bowtie2.yaml"
    container: "docker://ghcr.io/dartmouth-data-analytics-core/bowtie2:2.0"
    shell: "bowtie2-build {input} tests/fixtures/refs/tiny_genome"
