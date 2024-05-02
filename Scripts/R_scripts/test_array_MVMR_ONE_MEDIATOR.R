#!/cluster/software/EASYBUILD/R/4.1.0-foss-2021a/bin/Rscript
#SBATCH --job-name=MVMR_single_mediator
#SBATCH --output=./Slurm/Output/MVMR_single_mediator_%A_%a.out
#SBATCH --error=./Slurm/Error/MVMR_single_mediator_%A_%a.out
#SBATCH --account=p33
#SBATCH --time=1:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=50G
#SBATCH --array=1-5


source("Scripts/R_scripts/Lib/functions_for_logistics.R")
source("Scripts/R_scripts/Lib/functions_for_UVMR.R")
source("Scripts/R_scripts/Lib/functions_for_multivariate_MR_JB.R")
library(tidyverse)
library(TwoSampleMR)
library(MVMR)
library(parallel)


df_MDD <- read_tsv("./Data/Sumstats/MDD_Als_2023_sumstats") %>% 
  mutate(rsid = RSID, pval = P) %>% 
  mutate(phenotype = "MDD_Als_2023")

mediator_list <- get_MDD_var_list_without_CVD()
mediator <- mediator_list[1]
CVD <- get_MDD_var_list_only_updated_CVD()[1]
fit_MVMR_individual_mediator(df_MDD, mediator = mediator, CVD = CVD)


# 
# saveRDS(ivw, paste0("Results/Chunk_data/MVMR_single_mediators/MVMR_single_mediator_", CVD, "2.rds"))


