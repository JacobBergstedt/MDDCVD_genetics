#!/bin/bash
#SBATCH --job-name=mixer
#SBATCH --output=Slurm/Output/Mixer_plot_univariate.out
#SBATCH --error=Slurm/Error/Mixer_plot_univariate.out
#SBATCH --account=p33_tsd
#SBATCH --time=1:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=7G

module purge
module load singularity/3.7.3
export COMORMENT=/cluster/projects/p33/github/comorment
export SINGULARITY_BIND="$COMORMENT/mixer/reference:/REF:ro"
export SIF=$COMORMENT/mixer/singularity
export MIXER_COMMON_ARGS="--ld-file /REF/ldsc/1000G_EUR_Phase3_plink/1000G.EUR.QC.@.run4.ld --bim-file /REF/ldsc/1000G_EUR_Phase3_plink/1000G.EUR.QC.@.bim --threads 16"
export REP="rep${SLURM_ARRAY_TASK_ID}"
export EXTRACT="--extract /REF/ldsc/1000G_EUR_Phase3_plink/1000G.EUR.QC.prune_maf0p05_rand2M_r2p8.${REP}.snps"
export PYTHON="singularity exec --home ${PWD}:/home ${SIF}/mixer.sif python"

# AF CAD_2022 HF_2019_V2 PAD_2021_V2 Stroke_Gigastroke_2022 CRP_CHARGEUKB_2022 Educational_attainment LDL_GLGC_2021 PP Smoking DBP SBP HDL_GLGC_2021 TG_GLGC_2021 IL6 TC_GLGC_2021 TIID_DIAMANTE_2022 CM_2021_V2 BMI Sleep_duration_Doherty_2018_V2 Overall_physical_activity_Doherty_2018_V2 High_physical_activity_Wang_2022_V2 Loneliness NonHDL_GLGC_2021
for input in PAD
do

echo $input
$PYTHON /tools/mixer/precimed/mixer_figures.py combine --json Results/Mixer/MDD_Als_2023/${input}.test.rep@.json  --out Results/Mixer/MDD_Als_2023/Combined/${input}.test.combined
$PYTHON /tools/mixer/precimed/mixer_figures.py combine --json Results/Mixer/MDD_Als_2023/${input}.fit.rep@.json  --out Results/Mixer/MDD_Als_2023/Combined/${input}.fit.combined

$PYTHON /tools/mixer/precimed/mixer_figures.py one \
    --json Results/Mixer/MDD_Als_2023/Combined/${input}.test.combined.json \
    --out Results/Mixer/MDD_Als_2023/Output/${input}.test.output \
    --statistic mean std
    
$PYTHON /tools/mixer/precimed/mixer_figures.py one \
    --json Results/Mixer/MDD_Als_2023/Combined/${input}.fit.combined.json \
    --out Results/Mixer/MDD_Als_2023/Output/${input}.fit.output \
    --statistic mean std

done

# tar -zvcf /tsd/p33/data/durable/file-export/mixer_results_V2.5.tgz  Results/Mixer/Output/*AF* Results/Mixer/Output/*CAD_2022* Results/Mixer/Output/*Cardioembolic_stroke_Gigastroke_2022* Results/Mixer/Output/*Ischemic_stroke_Gigastroke_2022* Results/Mixer/Output/*Large_artery_stroke_Gigastroke_2022* Results/Mixer/Output/*Small_vessel_stroke_Gigastroke_2022* Results/Mixer/Output/*HF_2019_V2* Results/Mixer/Output/*PAD_2021_V2* Results/Mixer/Output/*Stroke_Gigastroke_2022*