library(tidyverse)
library(glue)
source("Scripts/R_scripts/Lib/functions_for_logistics.R")

res_individual_mediators <- glue("Results/Chunk_data/MVMR_single_mediators/MVMR_single_mediator_{get_MDD_var_list_only_updated_CVD()}_3.rds") %>% 
  map_dfr(readRDS) %>% 
  add_column(Analysis = "Individual")

res_individual_mediators$Covariate[!grepl("\\+", res_individual_mediators$model)] <- "None"


res_grouped_mediators <- glue("Results/Chunk_data/MVMR_grouped_mediators/MVMR_grouped_mediators_{get_MDD_var_list_only_updated_CVD()}3.rds") %>% 
  map_dfr(readRDS) %>% 
  add_column(Analysis = "Grouped")

res_only_MDD <- readRDS("Results/Chunk_data/MVMR_single_mediators/MVMR_single_only_MDD.rds") %>% 
  add_column(Analysis = "Only MDD")

res <- as_tibble(bind_rows(res_individual_mediators, res_grouped_mediators, res_only_MDD))

saveRDS(res, "Results/MVMR_results2.rds")
saveRDS(res, paste0(Sys.getenv("export"),"/MVMR_results2.rds"))

