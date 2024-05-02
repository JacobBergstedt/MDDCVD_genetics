library(tidyverse)
source("Scripts/R_scripts/Lib/functions_for_logistics.R")

path_filtered <- "/cluster/p/p33/cluster/users/johannap/ziyan/genomicSEMresults/CVD_MDD_GWAS_2_Q0.05.txt"
path <- "/cluster/p/p33/cluster/users/johannap/ziyan/genomicSEMresults/CVD_MDD_GWAS_2.txt"
keys <- get_keys()
MDDCVD_Q.0.05 <- read_delim(path_filtered) %>% 
  select(-i) %>% 
  rename(Z = Z_Estimate,
         P = Pval_Estimate,
         B = est) %>% 
  mutate(N = 17283 + 137371)

MDDCVD <- read_delim(path) %>% 
  select(-i) %>% 
  rename(Z = Z_Estimate,
         P = Pval_Estimate,
         B = est) %>% 
  mutate(N = 17283 + 137371)

write_tsv(MDDCVD, "Data/Sumstats/MDDCVD_latent_factor_sumstats")
write_tsv(MDDCVD_Q.0.05, "Data/Sumstats/MDDCVD_Q0.05_latent_factor_sumstats")