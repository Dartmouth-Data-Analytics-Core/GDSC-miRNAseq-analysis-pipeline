#!/bin/bash

#SBATCH --job-name=mirna
#SBATCH --nodes=1
#SBATCH --partition=preempt1
#SBATCH --account=dac
#SBATCH --mail-user=f007qps@dartmouth.edu
#SBATCH --time=60:00:00
#SBATCH --mail-type=FAIL
#SBATCH --output=mirna_%j.log

source /optnfs/common/miniconda3/etc/profile.d/conda.sh
conda activate /dartfs/rc/nosnapshots/G/GMBSR_refs/envs/snakemake

mkdir -p logs

snakemake -s \
    Snakefile \
    --configfile prebuilt_configs/zebrafish_config.yaml \
    --profile cluster_profile \
    -T 2 \
    --use-conda \
    --conda-frontend mamba \
    --conda-prefix /dartfs/rc/nosnapshots/G/GMBSR_refs/envs/miRNAseq
