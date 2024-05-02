library(tidyverse)
library(fs)

paths_uni <- fs::dir_ls("Results/Chunk_data/LAVA_V2.5/", glob = "*LAVA_MDD_uni*")
paths_biv <- fs::dir_ls("Results/Chunk_data/LAVA_V2.5/", glob = "*LAVA_MDD_biv*")
paths_pcor <- fs::dir_ls("Results/Chunk_data/LAVA_V2.5/", glob = "*LAVA_MDD_pcor*")

res_uni <- map_dfr(paths_uni, read_tsv)
res_biv <- map_dfr(paths_biv, read_tsv)
res_pcor <- map_dfr(paths_pcor, read_tsv)

saveRDS(res_uni, "/tsd/p33/data/durable/file-export/local_heritability_MDD_CVD_V2.5.rds")
saveRDS(res_biv, "/tsd/p33/data/durable/file-export/local_genetic_correlation_MDD_CVD_V2.5.rds")
saveRDS(res_pcor, "/tsd/p33/data/durable/file-export/local_genetic_correlation_MDD_CVD_adjusted_V2.5.rds")

saveRDS(res_uni, "Results/local_heritability_MDD_CVD_V2.5.rds")
saveRDS(res_biv, "Results/local_genetic_correlation_MDD_CVD_V2.5.rds")
saveRDS(res_pcor, "Results/local_genetic_correlation_MDD_CVD_adjusted_V2.5.rds")
