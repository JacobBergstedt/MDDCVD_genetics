library(tidyverse)
library(LAVA)
loci <- read.loci("Data/REF/LAVA/blocks_s2500_m25_f1_w200.GRCh37_hg19.locfile")
n_loc = nrow(loci)
split <- kronecker(1:100, rep(1, 25))[1:2495]
loci_list <- loci %>% split(f = split)
map2(loci_list, names(loci_list), ~ write.table(.x, paste0("./Data/TMP_data/LAVA_split/LAVA_loci_", .y, ".locfile") , quote = FALSE, row.names = FALSE))



loci <- read.loci("/cluster/projects/p33/github/comorment/containers/reference/")