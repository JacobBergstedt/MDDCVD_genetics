
library(TwoSampleMR)
library(tidyverse)
library(parallel)
source("./Scripts/R_scripts/Lib/functions_for_UVMR.R")
source("Scripts/R_scripts/Lib/functions_for_logistics.R")

# Params
MDD_trait <- "MDD_Als_2023_noUKB"
p_clumping <- 5e-8
r2_clumping <- 0.001


keys <- get_keys() %>% 
  mutate(units = if_else(Type == "Binary", "log odds", NA_character_)) %>% 
  mutate(prevalence = if_else(Type == "Binary", N_cases / N, NA_real_)) %>% 
  mutate(prevalence = if_else(grepl("MDD", Trait), 0.15, prevalence))

CVD_factor_set <- get_MDD_var_list_updated_CVD()


# MDD exposure ------------------------------------------------------------


MDD_input <- prepare_exposure(MDD_trait, keys = keys)

res_MDD_CVD_steiger <- mclapply(CVD_factor_set, 
                                function(outcome) UVMR(MDD_input, keys = keys, outcome = outcome, steiger_filter = TRUE, p_clumping = 5e-8, r2_clumping = 0.001), 
                                mc.cores = 12)

res_MDD_CVD_not_steiger <- mclapply(CVD_factor_set, 
                                    function(outcome) UVMR(MDD_input, keys = keys, outcome = outcome, steiger_filter = FALSE, p_clumping = 5e-8, r2_clumping = 0.001), 
                                    mc.cores = 12)


# CVD exposures -----------------------------------------------------------


CVD_inputs <- mclapply(CVD_factor_set, prepare_exposure, keys = keys, mc.cores = 24)

res_CVD_MDD_steiger <- mclapply(CVD_inputs, 
                                function(exposure_input) UVMR(exposure_input, keys = keys, outcome = MDD_trait, steiger_filter = TRUE, p_clumping = 5e-8, r2_clumping = 0.001), 
                                mc.cores = 24)

res_CVD_MDD_not_steiger <- mclapply(CVD_inputs, 
                                    function(exposure_input) UVMR(exposure_input, keys = keys, outcome = MDD_trait, steiger_filter = FALSE, p_clumping = 5e-8, r2_clumping = 0.001), 
                                    mc.cores = 24)


# Collect results ---------------------------------------------------------


res <- bind_rows(map_dfr(res_MDD_CVD_steiger, ~ .$res), 
                 map_dfr(res_MDD_CVD_not_steiger, ~ .$res),
                 map_dfr(res_CVD_MDD_steiger, ~ .$res),
                 map_dfr(res_CVD_MDD_not_steiger, ~ .$res)) %>% as_tibble()

pleiotropy <- bind_rows(add_column(map_dfr(res_MDD_CVD_steiger, ~ .$pleiotropy), Steiger_filtered = TRUE), 
                        add_column(map_dfr(res_MDD_CVD_not_steiger, ~ .$pleiotropy), Steiger_filtered = FALSE),
                        add_column(map_dfr(res_CVD_MDD_steiger, ~ .$pleiotropy), Steiger_filtered = TRUE), 
                        add_column(map_dfr(res_CVD_MDD_not_steiger, ~ .$pleiotropy), Steiger_filtered = FALSE)) %>% as_tibble()

heterogeneity <- bind_rows(add_column(map_dfr(res_MDD_CVD_steiger, ~ .$heterogeneity), Steiger_filtered = TRUE), 
                           add_column(map_dfr(res_MDD_CVD_not_steiger, ~ .$heterogeneity), Steiger_filtered = FALSE),
                           add_column(map_dfr(res_CVD_MDD_steiger, ~ .$heterogeneity), Steiger_filtered = TRUE), 
                           add_column(map_dfr(res_CVD_MDD_not_steiger, ~ .$heterogeneity), Steiger_filtered = FALSE)) %>% as_tibble()

# Write results -----------------------------------------------------------


write_tsv(res, "Results/MR/MDD_Als_2023/UVMR_leave_UKB_out_results.tsv")
write_tsv(pleiotropy, "Results/MR/MDD_Als_2023/UVMR_leave_UKB_out_pleiotropy.tsv")
write_tsv(heterogeneity, "Results/MR/MDD_Als_2023/UVMR_leave_UKB_out_heterogeneity.tsv")

write_tsv(res, paste0(Sys.getenv("export"), "/UVMR_results_Als_2023_leave_UKB_out.tsv"))
write_tsv(pleiotropy, paste0(Sys.getenv("export"), "/pleiotropy_Als_2023_leave_UKB_out.tsv"))
write_tsv(heterogeneity, paste0(Sys.getenv("export"), "/heterogeneity_Als_2023_leave_UKB_out.tsv"))
