#!/cluster/software/EASYBUILD/R/4.1.0-foss-2021a/bin/Rscript
library(tidyverse)
file <- fs::dir_ls("./Data/Munged_sumstats", glob = "*.sumstats")

trait_combs <- combn(file, m = 2, simplify = FALSE) %>% 
  map_dfr(~ tibble(T1 = .[[1]], T2 = .[[2]]))

write_delim(trait_combs, 
            "./Data/TMP_data/comorment_corr_combs.tsv", 
            col_names = FALSE)
