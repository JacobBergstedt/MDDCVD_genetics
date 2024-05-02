library(tidyverse)
source("Scripts/R_scripts/Lib/functions_for_logistics.R")
# traits <- get_MDD_var_list_only_updated_CVD()
traits <- "PAD"
paths <- paste0("Data/Sumstats/", traits, "_sumstats")

filter_mhc <- function(x) {
  
  if ("POS" %in% colnames(x)) {
    
    x <- x %>% rename(BP = POS)
    
  }
  
  if (!"Z" %in% colnames(x)) {
    
    x$Z <- x$B / x$SE
    
  }
  
  if ("RSID" %in% colnames(x)) {
    
    
    x <- rename(x, SNP = RSID)
    
    
  }
  
  if (!any(is.na(x$N_EFF))) x <- x %>% select(-N) %>% rename(N = N_EFF)
  
  x %>% 
    select(SNP, CHR, BP, A1, A2, N, Z ) %>% 
    filter(!is.na(SNP), !is.na(CHR), !is.na(BP), !is.na(A1), !is.na(A2), !is.na(N), !is.na(Z)) %>% 
    filter(!CHR %in% c("X", "Y", "M", "Mt"))
  

}


p <- map(paths, read_tsv) %>%
  map(filter_mhc) %>%
  map2(traits, ~ write_tsv(.x, file = paste0("Data/Sumstats_mixer_with_HLA/", .y, "_sumstats")))





