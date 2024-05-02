library(tidyverse)
library(fs)

paths_uni <- fs::dir_ls("Results/Chunk_data/LAVA_V2.5/", glob = "*LAVA_CM_uni*")
paths_biv <- fs::dir_ls("Results/Chunk_data/LAVA_V2.5/", glob = "*LAVA_CM_biv*")

res_uni <- map_dfr(paths_uni, read_tsv)
res_biv <- map_dfr(paths_biv, read_tsv)

p <- res_uni %>% filter(phen == "CM_2021_V2")
o <- res_uni %>% filter(phen != "CM_2021_V2")
res <- inner_join(o, p, by = c("locus", "chr", "start", "stop", "n.snps", "n.pcs")) %>% 
  rename(phen1 = phen.x, phen2 = phen.y, p.uni.1 = p.x, p.uni.2 = p.y, h2.uni.1 = h2.obs.x, h2.uni.2 = h2.obs.y) %>% 
  inner_join(res_biv)


saveRDS(res_uni, "/tsd/p33/data/durable/file-export/local_heritability_CM_MENT.rds")
saveRDS(res_biv, "/tsd/p33/data/durable/file-export/local_genetic_correlation_CM_MENT.rds")
saveRDS(res, "/tsd/p33/data/durable/file-export/local_genetic_correlation_CM_MENT_adjusted.rds")

saveRDS(res_uni, "Results/local_heritability_CM_MENT.rds")
saveRDS(res_biv, "Results/local_genetic_correlation_CM_MENT.rds")
