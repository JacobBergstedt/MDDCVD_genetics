#!/bin/bash
#SBATCH --job-name=SBayesS
#SBATCH --output=Slurm/Output/SBayesS_%A_%a.out
#SBATCH --error=Slurm/Error/SBayesS_%A_%a.out
#SBATCH --account=p33_norment
#SBATCH --cpus-per-task=5
#SBATCH --mem-per-cpu=10000M
#SBATCH --time=24:00:00
#SBATCH --array=1-24

trait_ID=${SLURM_ARRAY_TASK_ID}
trait=$(cat Data/TMP_data/keys_SBayesS.tsv | cut  -f2 | tail -n +2 | head -$trait_ID | tail -1)
path="Data/Sumstats_COJO/${trait}"
echo $trait
echo $path

echo $(date)

/cluster/p/p33/cluster/users/jacobb/gctb_2.05beta_Linux/gctb --sbayes S \
     --mldm Data/REF/ukbEURu_hm3_sparse/ukbEURu_hm3_shrunk_sparse/ukbEURu_hm3_sparse_mldm_list.txt \
     --pi 0.05 \
     --hsq 0.10 \
     --gwas-summary ${path} \
     --chain-length 15000 \
     --burn-in 5000 \
     --out Results/SBayesS/SBayesS${trait}

echo $(date)
