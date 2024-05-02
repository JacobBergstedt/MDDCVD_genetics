#!/bin/bash
#SBATCH --job-name=munge_sumstats
#SBATCH --output=./Slurm/Output/munge_sumstats.out
#SBATCH --error=./Slurm/Error/munge_sumstats.out
#SBATCH --account=p33_tsd
#SBATCH --time=5:00
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=10G

# module purge
# module load singularity/3.7.3
# export COMORMENT=/cluster/projects/p33/github/comorment
# export SIF=/cluster/projects/p33/github/comorment/ldsc/containers
# export SINGULARITY_BIND="/cluster/projects/p33/github/comorment/containers/reference:/REF:ro"

munge_fun () {
  local file=$1
  echo "TRAIT:""${file}"
  name=${file##*/}
  echo "${name}"
  singularity exec --home $PWD $COMORMENT/ldsc/containers/ldsc.sif python /tools/ldsc/munge_sumstats.py \
   --sumstats "./Data/Sumstats/""${name}" \
   --out "./Data/Munged_sumstats/""${name}""_munged" \
   --chunksize 100000 \
   --signed-sumstats B,0 \
   --merge-alleles "./Data/REF/w_hm3.snplist"
}

# export -f munge_fun
# parallel --jobs 10 munge_fun ::: ./Data/Sumstats/*_sumstats
#
# for sumstats in ./Data/Sumstats/MDDCVD_latent_factor_sumstats ./Data/Sumstats/MDDCVD_latent_factor_Q0.05_sumstats ./Data/Sumstats/MDDCVD_Als_2023_latent_factor_sumstats
# do
#   munge_fun ${sumstats}
# done


# for file in ./Data/Munged_sumstats/*_sumstats*gz
# do
#   new_file=${file::-3}
#   mv  "${file}" "${new_file}"
# done

# Process individual sumstat files
# munge_fun Data/Sumstats/ADHD_2019_sumstats
# munge_fun Data/Sumstats/Insomnia_2019_sumstats
# munge_fun Data/Sumstats/PTSD_2019_sumstats
# munge_fun Data/Sumstats/SCZ_2022_sumstats
# munge_fun Data/Sumstats/BIP_2021_sumstats
# munge_fun Data/Sumstats/AlcDep_2018_sumstats
# munge_fun Data/Sumstats/ANX_2020_sumstats
# munge_fun ./Data/Sumstats/MDD_Als_2023_sumstats

munge_fun ./Data/Sumstats/Sleep_duration_self_report_2019_sumstats

# munge_fun ./Data/Sumstats/Sleep_duration_short_self_report_2019_sumstats
# munge_fun ./Data/Sumstats/Sleep_duration_long_self_report_2019_sumstats

for file in ./Data/Munged_sumstats/Sleep_duration_self_report_2019_sumstats_munged.sumstats.gz 
do
  new_file=${file::-3}
  mv  "${file}" "${new_file}"
done