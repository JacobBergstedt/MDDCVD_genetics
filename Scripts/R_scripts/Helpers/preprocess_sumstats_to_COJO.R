library(tidyverse)
source("Scripts/R_scripts/Lib/functions_for_logistics.R")
keys <- get_keys() %>% 
  filter(Trait %in% get_MDD_var_list_updated_CVD())
keys_MDD <- get_keys() %>% filter(Trait == "MDD_Als_2023")
keys_PAD <- get_keys() %>% filter(Trait == "PAD_2021")

create_cojo_file <- function(x, i) {
  
  # HLA positions taken from https://www.ncbi.nlm.nih.gov/grc/human/regions/MHC?asm=GRCh37
  
  print(x$Trait)
  
  df <- read_tsv(x$Output_path) %>% 
    filter(!(CHR == 6 & (POS > 28477797 & POS < 33448354)))
  
  nrow(df)
  
  if (x$Trait %in% c("PP", "DBP", "SBP", "BMI", "CRP_CHARGEUKB_2022", "CRP_CHARGEUKB_2022", "CM_2021_V2")) {
    
    df <- df %>% select(SNP = RSID, A1, A2, freq = EAF_1KG, b = B, se = SE, p = P, N)
    
  } else {
    
    df <- df %>% select(SNP = RSID, A1, A2, freq = EAF, b = B, se = SE, p = P, N)
    
  }
  
  write_tsv(df, paste0("Data/Sumstats_COJO/", x$Trait))
  
}

keys %>%  
  rowwise() %>% 
  group_map(create_cojo_file)

create_cojo_file(keys_PAD, NULL)
