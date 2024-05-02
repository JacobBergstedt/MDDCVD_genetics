library(tidyverse)
library(glue)
source("Scripts/R_scripts/Lib/functions_for_logistics.R")
keys <- get_keys() %>% 
  filter(Trait %in% c("MDD_Als_2023", get_MDD_var_list_updated_CVD())) %>% 
  filter(Trait != "PAD_2021_V2")

res <- set_names(glue("Results/SBayesS/SBayesS{keys$Trait}.parRes"), keys$Trait) %>% 
  map(read.table) %>% 
  map_dfr(rownames_to_column, var = "Param", .id = "Trait")
  


write_tsv(res, "Results/SBayesS.tsv")
