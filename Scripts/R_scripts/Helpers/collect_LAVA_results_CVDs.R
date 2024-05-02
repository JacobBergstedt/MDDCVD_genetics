library(tidyverse)
library(fs)
library(glue)
source("Scripts/R_scripts/Lib/functions_for_logistics.R")

merge_uni_and_biv <- function(x, res_uni, res_biv) {
  
  p <- res_uni[[x]] %>% filter(phen == x)
  o <- res_uni[[x]] %>% filter(phen != x)
  inner_join(o, p, by = c("locus", "chr", "start", "stop", "n.snps", "n.pcs")) %>% 
    rename(phen1 = phen.x, phen2 = phen.y, p.uni.1 = p.x, p.uni.2 = p.y, h2.uni.1 = h2.obs.x, h2.uni.2 = h2.obs.y) %>% 
    inner_join(res_biv[[x]])
  
}


CVDs <- get_MDD_var_list_only_updated_CVD()
paths_uni <- map(set_names(CVDs, CVDs), ~ glue("Results/Chunk_data/LAVA_CVD/LAVA_uni_{.}_{1:100}.tsv"))
paths_biv <- map(set_names(CVDs, CVDs), ~ glue("Results/Chunk_data/LAVA_CVD/LAVA_biv_{.}_{1:100}.tsv"))

res_uni <- map(paths_uni, ~ map_dfr(., read_tsv))
res_biv <- map(paths_biv, ~ map_dfr(., read_tsv))

res <- map_dfr(CVDs, ~ merge_uni_and_biv(., res_uni, res_biv))

saveRDS(res, "Results/local_genetic_heritability_correlation_CVDs.rds")
saveRDS(res, paste0(Sys.getenv("export"), "/", "local_genetic_heritability_correlation_CVDs.rds"))
