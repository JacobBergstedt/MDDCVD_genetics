#!/cluster/software/EASYBUILD/R/4.1.0-foss-2021a/bin/Rscript
#SBATCH --job-name=MVMR_grouped_mediator
#SBATCH --output=./Slurm/Output/MVMR_grouped_mediator_%A_%a.out
#SBATCH --error=./Slurm/Error/MVMR_grouped_mediator_%A_%a.out
#SBATCH --account=p33
#SBATCH --time=1:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=25G
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

metabolic <- c("TG_GLGC_2021_sumstats",
               "HDL_GLGC_2021_sumstats",
               "NonHDL_GLGC_2021_sumstats",
               "TIID_DIAMANTE_2022_sumstats",
               "BMI_sumstats")

psychosocial <- c("High_physical_activity_Wang_2022_V2_sumstats",
                  "Loneliness_sumstats",
                  "Educational_attainment_sumstats",
                  "Smoking_sumstats")

inflammation <- c("CRP_CHARGEUKB_2022_sumstats")

cm <- c("CM_2021_V2_sumstats")

mediator_list <- list(metabolic = metabolic, psychosocial = psychosocial, inflammation = inflammation,cm = cm)
CVD <- get_MDD_var_list_only_updated_CVD()[1]

m <- fit_MVMR(df_MDD, mediators = psychosocial, CVD = CVD, covariate_group = "psychosocial")



# saveRDS(ivw, paste0("Results/Chunk_data/MVMR_grouped_mediators/MVMR_grouped_mediators_", CVD, "3.rds"))


