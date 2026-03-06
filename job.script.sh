#!/bin/bash

#SBATCH --job-name=mirna
#SBATCH --nodes=1
#SBATCH --partition=standard
#SBATCH --mail-user=f007qps@dartmouth.edu
#SBATCH --time=60:00:00
#SBATCH --mail-type=FAIL
#SBATCH --output=mirna_%j.log

source /optnfs/common/miniconda3/etc/profile.d/conda.sh
conda activate /dartfs/rc/nosnapshots/G/GMBSR_refs/envs/snakemake

mkdir -p slurm_logs

snakemake -s \
    Snakefile \
    --configfile prebuilt_configs/human_config.yaml \
    --profile cluster_profile \
    -T 2 \
    --use-conda \
    --conda-frontend conda \
    --conda-prefix /dartfs/rc/nosnapshots/G/GMBSR_refs/envs/miRNAseq
