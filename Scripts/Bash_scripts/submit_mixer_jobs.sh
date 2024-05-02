#!/bin/bash

# AF CAD_2022 HF_2019_V2 PAD_2021_V2 Stroke_Gigastroke_2022 
# CRP_CHARGEUKB_2022 Educational_attainment LDL_GLGC_2021 PP Smoking DBP SBP
# HDL_GLGC_2021 TG_GLGC_2021 IL6 TC_GLGC_2021 TIID_DIAMANTE_2022 CM_2021_V2
# BMI Sleep_duration_Doherty_2018_V2 Overall_physical_activity_Doherty_2018_V2 High_physical_activity_Wang_2022_V2 Loneliness NonHDL_GLGC_2021
# HF_2019_V2 PAD_2021_V2 Stroke_Gigastroke_2022 AF

for input in PAD
do

  sbatch Scripts/Bash_scripts/array_mixer_MDD_Als_2023_include_HLA.sh "${input}"

done

