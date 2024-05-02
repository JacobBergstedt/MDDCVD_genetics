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


CVD_list <- get_MDD_var_list_only_updated_CVD()
ivw <- mclapply(CVD_list, 
                function(CVD) fit_MVMR_only_MDD(df_MDD, keys = keys, CVD = CVD), 
                mc.cores = 4) %>% 
  bind_rows()

saveRDS(ivw, paste0("Results/Chunk_data/MVMR_single_mediators/MVMR_single_only_MDD.rds"))


