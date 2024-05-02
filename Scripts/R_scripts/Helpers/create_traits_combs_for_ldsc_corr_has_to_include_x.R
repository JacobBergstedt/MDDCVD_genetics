#!/cluster/software/EASYBUILD/R/4.1.0-foss-2021a/bin/Rscript
library(tidyverse)
file <- fs::dir_ls("./Data/Munged_sumstats", glob = "*.sumstats")

trait_combs <- combn(file, m = 2, simplify = FALSE) %>% 
  map_dfr(~ tibble(T1 = .[[1]], T2 = .[[2]]))

these <- c("Sleep_duration_short_self_report_2019", "Sleep_duration_long_self_report_2019")
include <- paste0("./Data/Munged_sumstats/", these, "_sumstats_munged.sumstats")

trait_combs <- trait_combs %>% 
  filter(T1 %in% include | T2 %in% include)

write_delim(trait_combs, 
            paste0("./Data/TMP_data/comorment_corr_combs_has_to_include_sleep_duration_traits.tsv"), 
            col_names = FALSE)