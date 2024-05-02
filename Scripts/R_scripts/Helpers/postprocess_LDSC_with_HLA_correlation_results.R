library(tidyverse)
library(vroom)
library(corrr)
library(fs)
source("Scripts/R_scripts/Lib/helpers.R")

p <- map_dfr(dir_ls("./Results/Chunk_data/LDSC_corr_including_HLA/", glob = "*sumstats"), read_table) %>%
  filter(grepl("MDD_Als_2023", p1) | grepl("MDD_Als_2023", p2)) %>% 
  mutate(Trait1 =  str_extract(p1, "(?<=Munged_sumstats/).*(?=_sumstats)"),
         Trait2 =  str_extract(p2, "(?<=Munged_sumstats/).*(?=_sumstats)"))

write_tsv(p, "Results/LDSC_corr/ldsc_correlations_with_HLA_Als_2023.tsv")
write_tsv(p, "/tsd/p33/data/durable/file-export/ldsc_correlations_with_HLA_Als_2023.tsv")