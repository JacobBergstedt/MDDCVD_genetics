#!/bin/bash
#SBATCH --job-name=V2_latentCVD
#SBATCH --account=p33
#SBATCH --time=1-12:00:00
#SBATCH --cpus-per-task=32
#SBATCH --mem-per-cpu=8000M

module load singularity/3.7.1

singularity exec --home $PWD $SIF/r.sif Rscript Scripts/R_scripts/latentCVD_GWAS_V2.5.R 22
