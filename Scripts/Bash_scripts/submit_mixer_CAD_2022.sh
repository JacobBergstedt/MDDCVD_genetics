#!/bin/bash
#SBATCH --job-name=mixer
#SBATCH --output=Slurm/Output/Mixer_CAD_2022_%A_%a.out
#SBATCH --error=Slurm/Error/Mixer_CAD_2022_%A_%a.out
#SBATCH --account=p33_norment
#SBATCH --time=24:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem-per-cpu=8000M
#SBATCH --array=1-20

module purge
module load singularity/3.7.3
export COMORMENT=/cluster/projects/p33/github/comorment
export SINGULARITY_BIND="$COMORMENT/mixer/reference:/REF:ro"
export SIF=$COMORMENT/mixer/singularity
export MIXER_COMMON_ARGS="--ld-file /REF/ldsc/1000G_EUR_Phase3_plink/1000G.EUR.QC.@.run4.ld --bim-file /REF/ldsc/1000G_EUR_Phase3_plink/1000G.EUR.QC.@.bim --threads 8"
export REP="rep${SLURM_ARRAY_TASK_ID}"
export EXTRACT="--extract /REF/ldsc/1000G_EUR_Phase3_plink/1000G.EUR.QC.prune_maf0p05_rand2M_r2p8.${REP}.snps"
export PYTHON="singularity exec --home ${PWD}:/home ${SIF}/mixer.sif python"

CAD_path="Data/Sumstats_without_HLA/CAD_2022_sumstats"
CAD_fit_out="Results/Mixer/CAD_2022/CAD_2022.fit.${REP}"
CAD_param="${CAD_fit_out}.json"
CAD_test_out="Results/Mixer/CAD_2022/CAD_2022.test.${REP}"

${PYTHON} /tools/mixer/precimed/mixer.py fit1 ${MIXER_COMMON_ARGS} ${EXTRACT} --trait1-file "${CAD_path}" --out ${CAD_fit_out}
${PYTHON} /tools/mixer/precimed/mixer.py test1 ${MIXER_COMMON_ARGS} --trait1-file "${CAD_path}" --load-params "${CAD_param}" --out "${CAD_test_out}"