#!/bin/bash
for trait in AF CAD_2022 HF_2019_V2 PAD_2021_V2 Stroke_Gigastroke_2022
do
  sbatch Scripts/R_scripts/array_local_genetic_correlation_CVD_trait_2023.R --target="${trait}"
done
