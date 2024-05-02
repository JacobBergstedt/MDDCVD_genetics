library(tidyverse)
library(fs)
library(parallel)

trait <- "MDDCVD_Q0.05_Als_2023_latent_factor.FUSION"
paths <- dir_ls("Results/TWAS", regexp = trait)
res1 <- mclapply(paths, read_tsv, mc.cores = 10) %>% 
  bind_rows() %>% 
  mutate(Trait = "MDDCVD_Q0.05_Als_2023_latent_factor")


trait <- "CVD_latent_factor.FUSION"
paths <- dir_ls("Results/TWAS", regexp = trait)
res2 <- mclapply(paths, read_tsv, mc.cores = 10) %>% 
  bind_rows() %>% 
  mutate(Trait = "CVD_latent_factor")

trait <- "MDD_Als_2023.FUSION"
paths <- dir_ls("Results/TWAS", regexp = trait)
res3 <- mclapply(paths, read_tsv, mc.cores = 10) %>% 
  bind_rows() %>% 
  mutate(Trait = "MDD_Als_2023")


trait <- "MDDCVD_Als_2023_latent_factor.FUSION"
paths <- dir_ls("Results/TWAS", regexp = trait)
res4 <- mclapply(paths, read_tsv, mc.cores = 10) %>% 
  bind_rows() %>% 
  mutate(Trait = "MDDCVD_Als_2023_latent_factor")


res <- bind_rows(res1, res2, res3, res4)
