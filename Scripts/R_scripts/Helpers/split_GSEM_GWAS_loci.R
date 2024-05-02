library(tidyverse)
MDDCVDsumstats <- readRDS("Data/Sumstats_for_GSEM_GWAS/EFA_CVD_MDD_V2.5_sumstats.rds")
nr_splits <- 300
idx <- rep_len(1:nr_splits, length.out = nrow(MDDCVDsumstats)) %>% sort()
MDDCVDsumstats$IDX <- idx

for (i in 1:nr_splits) {
  
  MDDCVDsumstats_sub <- MDDCVDsumstats %>% filter(IDX == i) %>% select(-IDX)
  saveRDS(MDDCVDsumstats_sub, paste0("Data/TMP_data/GSEM_split/MDD_CVD_GSEM_LOCI_SPLIT", i, ".rds"))
  
}

CVD_sumstats <- readRDS("Data/Sumstats_for_GSEM_GWAS/EFA_CVD_V2.5_sumstats.rds")
idx <- rep_len(1:nr_splits, length.out = nrow(CVD_sumstats)) %>% sort()
CVD_sumstats$IDX <- idx

for (i in 1:nr_splits) {
  
  CVD_sumstats_sub <- CVD_sumstats %>% filter(IDX == i) %>% select(-IDX)
  saveRDS(CVD_sumstats_sub, paste0("Data/TMP_data/GSEM_split/CVD_GSEM_LOCI_SPLIT", i, ".rds"))
  
}
