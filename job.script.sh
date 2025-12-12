#!/bin/bash

# Name of the job
#SBATCH --job-name=mirna

# Number of compute nodes
#SBATCH --nodes=1

# partition
#SBATCH --partition=standard

# account
#SBATCH --account=nccc

# Walltime (job duration)
#SBATCH --time=60:00:00

# Email notifications (comma-separated options: BEGIN,END,FAIL)
#SBATCH --mail-type=FAIL

source /optnfs/common/miniconda3/etc/profile.d/conda.sh

conda activate /dartfs/rc/nosnapshots/G/GMBSR_refs/envs/snakemake

snakemake -s Snakefile --profile cluster_profile -T 2 --use-conda --conda-frontend conda --conda-prefix /dartfs/rc/nosnapshots/G/GMBSR_refs/envs/miRNAseq
