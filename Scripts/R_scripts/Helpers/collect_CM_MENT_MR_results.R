library(fs)
library(tidyverse)

res <- map_dfr(fs::dir_ls("Results/MR/CM", glob = "UVMR_V2*main*"), read_tsv, .id = "Path")
write_tsv(res, "/tsd/p33/data/durable/file-export/UVMR_CM_MENT_V2.tsv")
