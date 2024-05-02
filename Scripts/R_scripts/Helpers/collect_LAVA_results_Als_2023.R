library(tidyverse)
library(fs)

paths_uni <- fs::dir_ls("Results/Chunk_data/LAVA_Als_2023/", glob = "*LAVA_uni_*")
paths_biv <- fs::dir_ls("Results/Chunk_data/LAVA_Als_2023/", glob = "*LAVA_biv*")

res_uni <- map_dfr(paths_uni, read_tsv)
res_biv <- map_dfr(paths_biv, read_tsv)

p <- res_uni %>% filter(phen == "MDD_Als_2023")
o <- res_uni %>% filter(phen != "MDD_Als_2023")
res <- inner_join(o, p, by = c("locus", "chr", "start", "stop", "n.snps", "n.pcs")) %>% 
  rename(phen1 = phen.x, phen2 = phen.y, p.uni.1 = p.x, p.uni.2 = p.y, h2.uni.1 = h2.obs.x, h2.uni.2 = h2.obs.y) %>% 
  inner_join(res_biv)

saveRDS(res_uni, "Results/local_heritability_Als_2023.rds")
saveRDS(res_biv, "Results/local_genetic_correlation_Als_2023.rds")
saveRDS(res, "Results/local_genetic_heritability_correlation_Als_2023.rds")

saveRDS(res_uni, paste0(Sys.getenv("export"), "/", "local_genetic_heritability_Als_2023.rds"))
saveRDS(res_biv, paste0(Sys.getenv("export"), "/", "local_genetic_correlation_Als_2023.rds"))
saveRDS(res, paste0(Sys.getenv("export"), "/", "local_genetic_heritability_correlation_Als_2023.rds"))