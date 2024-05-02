library(TwoSampleMR)
library(tidyverse)
library(parallel)
source("./Scripts/R_scripts/Lib/functions_for_UVMR.R")
source("Scripts/R_scripts/Lib/functions_for_logistics.R")

run_UVMR_for_CVDs <- function(input, CVD_factor_set) {
  
  out <- map(CVD_factor_set, ~ UVMR(input, outcome = ., r2_clumping = 0.1))
  names(out) <- CVD_factor_set
  out
}

CVD_factor_set <- get_MDD_var_list_only_updated_CVD()
exposures <- mclapply(set_names(get_ment_list(), get_ment_list()), prepare_exposure, mc.cores = 10)
res <- mclapply(exposures, function(input) run_UVMR_for_CVDs(input, CVD_factor_set), mc.cores = 10)



res_out <- map_dfr(res, function(x) map_dfr(x, function(y) y$res)) %>% as_tibble()
pleiotropy <- map_dfr(res, function(x) map_dfr(x, function(y) y$pleiotropy)) %>% as_tibble()
heterogeneity <- map_dfr(res, function(x) map_dfr(x, function(y) y$heterogeneity)) %>% as_tibble()


write_tsv(res_out, paste0(Sys.getenv("export"), "/UVMR_results_ment_CVD.tsv"))
write_tsv(pleiotropy, paste0(Sys.getenv("export"), "/pleiotropy_ment_CVD.tsv"))
write_tsv(heterogeneity, paste0(Sys.getenv("export"), "/heterogeneity_ment_CVD.tsv"))

