library(tidyverse)

# MDD no UKB --------------------------------------------------------------
mdd_noUKB_orig <- read_tsv("../processed_sumstats_ForJacob/PGC_MDD_2022_no_23andME_noUKB/raw/daner_pgc_mdd_ex.23aME.loo.no.UKBB.gz") %>%
  select(CHR, RSID = SNP, POS = BP, SE)

mdd_noUKB <- read_tsv("Data/Sumstats/MDD_noUKB_sumstats") %>% 
  left_join(mdd_noUKB_orig)
