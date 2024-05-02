library(tidyverse)
library(fs)
paths <- dir_ls("Results/Chunk_data/GSEM_CVD/")

error_chr <- function(x) {
  
  x$error <- as.character(x$error)
  x
  
}
est_sample_size <- function(x) {
  
  x <- x %>% filter(EAF <= 0.4, EAF >= 0.1)
  EAF <- x$EAF
  var <- 2 * EAF * (1 - EAF)
  mean(1 / x$SE ^ 2 * 1 / var)
  
}

df <- mclapply(paths, readRDS, mc.cores = 12) %>% 
  bind_rows() %>% 
  filter(!is.na(se_c)) %>% 
  select(RSID = SNP, CHR, POS = BP, EAF = MAF, A1, A2, B = est, SE = se_c, Z = Z_Estimate, P = Pval_Estimate, Q, Q_df, Q_pval) %>%
  mutate(N = est_sample_size(.))

write_tsv(df, "Data/Sumstats/CVD_latent_factor_sumstats")