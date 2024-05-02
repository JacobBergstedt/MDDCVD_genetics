#!/cluster/software/EASYBUILD/R/4.1.0-foss-2021a/bin/Rscript
#SBATCH --job-name=MVMR_grouped_mediator
#SBATCH --output=./Slurm/Output/MVMR_grouped_mediator_%A_%a.out
#SBATCH --error=./Slurm/Error/MVMR_grouped_mediator_%A_%a.out
#SBATCH --account=p33_norment
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

keys <- get_keys() %>% 
  mutate(units = if_else(Type == "Binary", "log odds", NA_character_)) %>% 
  mutate(prevalence = if_else(Type == "Binary", N_cases / N, NA_real_)) %>% 
  mutate(prevalence = if_else(grepl("MDD", Trait), 0.15, prevalence))

df_MDD <- read_tsv("./Data/Sumstats/MDD_Als_2023_sumstats") %>% 
  mutate(rsid = RSID, prevalence = 0.15, pval = P, phenotype = "MDD_Als_2023", units = "log odds")


metabolic <- c("TG_GLGC_2021",
               "HDL_GLGC_2021",
               "NonHDL_GLGC_2021",
               "TIID_DIAMANTE_2022",
               "BMI")

psychosocial <- c("High_physical_activity_Wang_2022_V2",
                  "Loneliness",
                  "Educational_attainment",
                  "Smoking")

inflammation <- "CRP_CHARGEUKB_2022"
cm <- "CM_2021_V2"

mediator_list <- list(metabolic = metabolic, psychosocial = psychosocial, inflammation = inflammation, cm = cm)
CVD <- get_MDD_var_list_only_updated_CVD()[as.numeric(Sys.getenv("SLURM_ARRAY_TASK_ID"))]


ivw <- mclapply(names(mediator_list), 
                function(str) fit_MVMR(df_MDD, keys = keys, mediators = mediator_list[[str]], CVD = CVD, covariate_group = str), 
                mc.cores = 4) %>% 
  bind_rows()

saveRDS(ivw, paste0("Results/Chunk_data/MVMR_grouped_mediators/MVMR_grouped_mediators_", CVD, "3.rds"))


