library(TwoSampleMR)
library(tidyverse)
library(parallel)
source("./Scripts/R_scripts/Lib/functions_for_UVMR.R")
source("Scripts/R_scripts/Lib/functions_for_logistics.R")

CVD_factor_set <- c(get_MDD_Als(), get_MDD_var_list_updated_CVD())
MDD_input <- prepare_exposure("MDDCVD_Als_2023_latent_factor_Q0.05")
res_MDD_CVD <- mclapply(CVD_factor_set, 
                        function(outcome) UVMR(MDD_input, outcome = outcome, r2_clumping = 0.001), 
                        mc.cores = 25)

CVD_inputs <- mclapply(c(get_MDD_Als(), get_MDD_var_list_without_CVD()), prepare_exposure, mc.cores = 20)
res_CVD_MDD <- mclapply(CVD_inputs, 
                        function(exposure_input) UVMR(exposure = exposure_input, outcome = "MDDCVD_Als_2023_latent_factor_Q0.05", r2_clumping = 0.001), 
                        mc.cores = 20)

res <- bind_rows(map_dfr(res_MDD_CVD, ~ .$res), map_dfr(res_CVD_MDD, ~ .$res)) %>% as_tibble()
pleiotropy <- bind_rows(map_dfr(res_MDD_CVD, ~ .$pleiotropy), map_dfr(res_CVD_MDD, ~ .$pleiotropy)) %>% as_tibble()
heterogeneity <- bind_rows(map_dfr(res_MDD_CVD, ~ .$heterogeneity), map_dfr(res_CVD_MDD, ~ .$heterogeneity)) %>% as_tibble()

write_tsv(res, "Results/MR/MDD_Als_2023/UVMR_latent.tsv")
write_tsv(pleiotropy, "Results/MR/MDD_Als_2023/UVMR_latent_pleiotropy.tsv")
write_tsv(heterogeneity, "Results/MR/MDD_Als_2023/UVMR_latent_heterogeneity.tsv")

write_tsv(res, paste0(Sys.getenv("export"), "/UVMR_results_Als_2023_latent.tsv"))
write_tsv(pleiotropy, paste0(Sys.getenv("export"), "/pleiotropy_Als_2023_latent.tsv"))
write_tsv(heterogeneity, paste0(Sys.getenv("export"), "/heterogeneity_Als_2023_latent.tsv"))