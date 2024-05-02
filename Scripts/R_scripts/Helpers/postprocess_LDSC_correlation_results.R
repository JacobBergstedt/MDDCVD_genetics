library(tidyverse)
library(vroom)
library(corrr)
library(fs)
source("Scripts/R_scripts/Lib/helpers.R")

p <- map_dfr(dir_ls("./Results/Chunk_data/LDSC_corr", glob = "*sumstats"), read_table) %>%
  filter(!grepl("MDDHowardCVD_latent_factor_munged.sumstats", p1), 
         !grepl("MDDHowardCVD_latent_factor_munged.sumstats", p2),
         !grepl("MDD_Howard_sumstats_munged.sumstats", p1),
         !grepl("MDD_Howard_sumstats_munged.sumstats", p2)) %>% 
  mutate(Trait1 =  str_extract(p1, "(?<=Munged_sumstats/).*(?=_sumstats)"),
         Trait2 =  str_extract(p2, "(?<=Munged_sumstats/).*(?=_sumstats)"))

sample_overlap_mat <- cor_mat_from_long_format(p, "gcov_int")
write.table(as.data.frame(sample_overlap_mat), "./Data/LDSC_sample_overlap_matrix_Als_2023.txt", quote = FALSE)
write_tsv(p, "Results/LDSC_corr/ldsc_correlations_Als_2023.tsv")
write_tsv(p, "/tsd/p33/data/durable/file-export/ldsc_correlations_Als_2023.tsv")