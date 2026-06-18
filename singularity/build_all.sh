#!/usr/bin/env bash
# Build all Singularity images from definition files.
# Run from singularity/defs/:
#  cd singularity/defs && bash ../build_all.sh

singularity build ../bbmap.sif bbmap.def
singularity build ../cutadapt.sif cutadapt.def
singularity build ../fastqc.sif fastqc.def
singularity build ../featurecounts.sif featurecounts.def
singularity build ../mirtop.sif mirtop.def
singularity build ../multiqc.sif multiqc.def
singularity build ../pcaplot.sif pcaplot.def
singularity build ../picard.sif picard.def
singularity build ../r_env.sif r_env.def
singularity build ../seqcluster.sif seqcluster.def
