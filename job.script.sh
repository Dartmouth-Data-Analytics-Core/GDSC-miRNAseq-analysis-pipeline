#!/bin/bash

#SBATCH --job-name=mirna                          
#SBATCH --nodes=1
#SBATCH --partition=standard
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=16  
#SBATCH --time=60:00:00
#SBATCH --mail-user=f007qps@dartmouth.edu
#SBATCH --mail-type=FAIL
#SBATCH --output=%x_%j.log
#========================================================#

#----- Specify Config (one of "human", "mouse", or "zebrafish", case-sensitive and needs to be in quotes.)
CONFIG="human"

#----- Environment information
CONDA_BASE="/optnfs/common/miniconda3"
SNAKEMAKE_ENV="/dartfs/rc/nosnapshots/G/GMBSR_refs/envs/snakemake"
CONDA_PREFIX_PATH="/dartfs/rc/nosnapshots/G/GMBSR_refs/envs/GDSC-Clover-Seq"
source "${CONDA_BASE}/etc/profile.d/conda.sh"
conda activate "${SNAKEMAKE_ENV}"

#----- LOGGER
echo "#------------------------ Initialization ------------------------#"
echo -e "Running GDSC-miRNASeq Pipeline for "${CONFIG}" with Snakemake $(snakemake --version) \n" 
echo "Job:        $SLURM_JOB_NAME"
echo "Job ID:     $SLURM_JOB_ID"
echo "Node:       $(hostname)"
echo "Start time: $(date)"
echo "Work dir:   $(pwd)"
echo "Conda base: $CONDA_BASE"
echo "Snakemake environment: $SNAKEMAKE_ENV"
echo "Using snakemake from: $(which snakemake)"
echo "Conda Prefix: $CONDA_PREFIX_PATH"
echo -e "#-------------------------------------------------------------# \n"
echo -e "SNAKEMAKE LOG:\n"

#----- Make slurm logs
mkdir -p slurm_logs

#----- Invoke Snakemake
snakemake -s \
    Snakefile \
    --configfile prebuilt_configs/"${CONFIG}"_config.yaml \
    --profile cluster_profile \
    -T 2 \
    --use-conda \
    --conda-frontend conda \
    --conda-prefix /dartfs/rc/nosnapshots/G/GMBSR_refs/envs/miRNAseq

#----- END
echo "End time: $(date)"
