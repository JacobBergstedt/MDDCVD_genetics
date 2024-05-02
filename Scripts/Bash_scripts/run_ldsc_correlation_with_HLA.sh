#!/bin/bash
#SBATCH --job-name=ldsc_correlation_include_HLA
#SBATCH --output=./Slurm/Output/ldsc_include_HLA_corr_out
#SBATCH --error=./Slurm/Error/ldsc_include_HLA_corr_out
#SBATCH --account=p33
#SBATCH --time=5:00:00
#SBATCH --cpus-per-task=2
#SBATCH --mem-per-cpu=10G

module purge
module load singularity/3.7.3

export COMORMENT=/cluster/projects/p33/github/comorment
export SIF=/cluster/projects/p33/github/comorment/ldsc/containers
export SINGULARITY_BIND="/cluster/projects/p33/github/comorment/containers/reference:/REF:ro"

for i in {1..87}
do
  p=$(head -n $i ./Data/TMP_data/comorment_corr_combs_has_to_include_MDD_Als_2023.tsv | tail -1)
  file1=$(echo "$p" | awk '{print $1}')
  file2=$(echo "$p" | awk '{print $2}')
  trait1=${file1##*/}
  trait2=${file2##*/}
  outpath="Results/Chunk_data/LDSC_corr_including_HLA/LDSC_corr_inc_HLA_""${trait1}""_""${trait2}"

  echo $file1
  echo $file2

  singularity exec --home $PWD:/home ${SIF}/ldsc.sif python /tools/ldsc/ldsc.py \
    --rg "${file1}","${file2}" \
    --ref-ld-chr Data/REF/eur_w_ld_chr_include_HLA/ \
    --w-ld-chr Data/REF/eur_w_ld_chr_include_HLA/ \
    --out "${outpath}"

  awk '/Summary/{flag=1; next} /Analysis/{flag=0} flag' "${outpath}"".log" > $outpath
done


# file1="./Data/Munged_sumstats/Loneliness_sumstats_munged.sumstats"
# file2="./Data/Munged_sumstats/MDD_Als_2023_sumstats_munged.sumstats"
# trait1=${file1##*/}
# trait2=${file2##*/}
# outpath="Results/Chunk_data/LDSC_corr_including_HLA/LDSC_including_HLA_corr_""${trait1}""_""${trait2}"
# 
# echo $file1
# echo $file2
# 
# singularity exec --home $PWD:/home ${SIF}/ldsc.sif python /tools/ldsc/ldsc.py \
#   --rg "${file1}","${file2}" \
#   --ref-ld-chr Data/REF/eur_w_ld_chr_include_HLA/ \
#   --w-ld-chr Data/REF/eur_w_ld_chr_include_HLA/ \
#   --out "${outpath}"
# 
# awk '/Summary/{flag=1; next} /Analysis/{flag=0} flag' "${outpath}"".log" > $outpath


# p=$(head -n 10 "Data/TMP_data/comorment_corr_combs.tsv" | tail -1)
# file1=$(echo "$p" | awk '{print $1}')
# file2=$(echo "$p" | awk '{print $2}')
# trait1=${file1##*/}
# trait2=${file2##*/}
# outpath="Data/TMP_data/LDSC_corr_""${trait1}""_""${trait2}"
# 
# echo $file1
# echo $file2
# 
# singularity exec --home $PWD:/home ${SIF}/ldsc.sif python /tools/ldsc/ldsc.py \
#   --rg "${file1}","${file2}" \
#   --ref-ld-chr Data/REF/eur_w_ld_chr_include_HLA/ \
#   --w-ld-chr Data/REF/eur_w_ld_chr_include_HLA/ \
#   --out "${outpath}"
# 
# awk '/Summary/{flag=1; next} /Analysis/{flag=0} flag' "${outpath}"".log" > $outpath
# 
