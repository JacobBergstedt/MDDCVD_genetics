#!/cluster/software/EASYBUILD/R/4.1.0-foss-2021a/bin/Rscript
#SBATCH --job-name=MVMR_single_mediator
#SBATCH --output=./Slurm/Output/MVMR_single_mediator_%A_%a.out
#SBATCH --error=./Slurm/Error/MVMR_single_mediator_%A_%a.out
#SBATCH --account=p33_norment
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

keys <- get_keys() %>% 
  mutate(units = if_else(Type == "Binary", "log odds", NA_character_)) %>% 
  mutate(prevalence = if_else(Type == "Binary", N_cases / N, NA_real_)) %>% 
  mutate(prevalence = if_else(grepl("MDD", Trait), 0.15, prevalence))


df_MDD <- read_tsv("./Data/Sumstats/MDD_Als_2023_sumstats") %>% 
  mutate(rsid = RSID, prevalence = 0.15, pval = P, phenotype = "MDD_Als_2023", units = "log odds")

mediator_list <- get_MDD_var_list_without_CVD()
CVD <- get_MDD_var_list_only_updated_CVD()[as.numeric(Sys.getenv("SLURM_ARRAY_TASK_ID"))]

print(CVD)
print(mediator_list)

ivw <- mclapply(mediator_list, 
                function(mediator) fit_MVMR_individual_mediator(df_MDD, keys = keys, mediator = mediator, CVD = CVD), 
                mc.cores = 4) %>% 
  bind_rows()

saveRDS(ivw,
        paste0("Results/Chunk_data/MVMR_single_mediators/MVMR_single_mediator_", CVD, "_3.rds"))


