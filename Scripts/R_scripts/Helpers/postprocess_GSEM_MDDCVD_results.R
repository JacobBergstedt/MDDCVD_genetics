library(tidyverse)
library(fs)
library(parallel)
paths_common <- dir_ls("Results/Chunk_data/GSEM_MDDCVD/Als_2023/",  regexp = "COMMON")
paths_independent <- dir_ls("Results/Chunk_data/GSEM_MDDCVD/Als_2023/",  regexp = "INDEPENDENT")

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

df <- mclapply(paths_common, readRDS, mc.cores = 12) %>% 
  flatten() %>% 
  map(error_chr) %>% 
  bind_rows() %>% 
  filter(lhs == "CVDMDD", op == "~", rhs == "SNP") %>% 
  as_tibble()

df_independent_pathway <- mclapply(paths_independent, readRDS, mc.cores = 12) %>% 
  flatten() %>% 
  map(error_chr) %>% 
  bind_rows()

df_independent_pathway_chisq <- select(df_independent_pathway, SNP, CHR, chisq_independent = chisq, chisq_pval_independent = chisq_pval, chisq_df_independent = chisq_df) %>% 
  distinct(SNP, CHR, .keep_all = TRUE)

df <- df %>% 
  left_join(df_independent_pathway_chisq) %>% 
  mutate(Q = chisq - chisq_independent, Q_df = chisq_df - chisq_df_independent) %>% 
  mutate(Q_pval = pchisq(Q, Q_df, lower.tail = FALSE)) %>% 
  select(RSID = SNP, CHR, POS = BP, EAF = MAF, A1, A2, B = est, SE, Z = Z_Estimate, P = Pval_Estimate, Q, Q_pval, Q_df) %>%
  mutate(N = est_sample_size(.))

df_filter_Q <- df %>% 
  filter(Q_pval > 0.05) %>% 
  mutate(N = est_sample_size(.))

write_tsv(df, "Data/Sumstats/MDDCVD_Als_2023_latent_factor_sumstats")
write_tsv(df_filter_Q, "Data/Sumstats/MDDCVD_Als_2023_latent_factor_Q0.05_sumstats")

