paths <- fs::dir_ls(path = "Results/Chunk_data/Lhc_MR/")
res <- map_dfr(paths, read_tsv, .id = "Path") %>% 
  mutate(X = "MDD", Y = str_extract(Path, "(?=MDD).*(?>_sumstats)") %>% str_sub(start = 4, end = -10))
saveRDS(res, "Results/res_lhc_MR.rds")
