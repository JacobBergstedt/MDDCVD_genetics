
library(tidyverse)

##################################################################################################
### Prep the clean sumstats files for the GWASs
### The software doesn't like having 2 effect columns, so remove Z-column

CAD <- read_tsv("Data/Sumstats/CAD_2022_sumstats") %>% 
  select(-OR) %>% 
  rename(MAF = EAF) %>%  
  write_tsv("Data/Sumstats_for_GSEM_GWAS/CAD_2022_sumstats.gsem")

MDD <- read_tsv("Data/Sumstats/MDD_PGC3_V2_sumstats") %>% 
  filter(!is.na(POS)) %>% 
  rename(MAF = EAF) %>%  
  write_tsv("Data/Sumstats_for_GSEM_GWAS/MDD_PGC3_V2_sumstats.gsem")

STR <- read_tsv("Data/Sumstats/Stroke_Gigastroke_2022_sumstats") %>% 
  select(-Z, -OR) %>% 
  rename(MAF = EAF) %>% 
  write_tsv("Data/Sumstats_for_GSEM_GWAS/Stroke_Gigastroke_2022_sumstats.gsem")

PAD <- read_tsv("Data/Sumstats/PAD_2021_V2_sumstats") %>% 
  select(-Z, -OR) %>% 
  rename(MAF = EAF) %>% 
  write_tsv("Data/Sumstats_for_GSEM_GWAS/PAD_2021_V2_sumstats.gsem")

HF <- read_tsv("Data/Sumstats/HF_2019_V2_sumstats") %>%
  select(-Z) %>% 
  rename(MAF = EAF) %>% 
  write_tsv("Data/Sumstats_for_GSEM_GWAS/HF_2019_V2_sumstats.gsem")


MDDCVD_files <- c("Data/Sumstats_for_GSEM_GWAS/CAD_2022_sumstats.gsem",
                  "Data/Sumstats_for_GSEM_GWAS/MDD_PGC3_V2_sumstats.gsem",
                  "Data/Sumstats_for_GSEM_GWAS/Stroke_Gigastroke_2022_sumstats.gsem",
                  "Data/Sumstats_for_GSEM_GWAS/PAD_2021_V2_sumstats.gsem",
                  "Data/Sumstats_for_GSEM_GWAS/HF_2019_V2_sumstats.gsem")

ref= "Data/REF/reference.1000G.maf.0.005.txt"
MDDCVD_names=c("CAD_2022", 
              "MDD_PGC3_V2",
              "Stroke_Gigastroke_2022",
              "PAD_2021_V2",
              "HF_2019_V2")

se.logit <- c(TRUE, TRUE, TRUE, TRUE, TRUE)

MDDCVD_sumstats <- sumstats(files = MDDCVD_files,
                            ref = ref,
                            trait.names = MDDCVD_names,
                            se.logit = se.logit)

saveRDS(MDDCVD_sumstats, "Data/Sumstats_for_GSEM_GWAS/EFA_CVD_MDD_V2.5_sumstats.rds")

CVD_files <- c("Data/Sumstats_for_GSEM_GWAS/CAD_2022_sumstats.gsem",
               "Data/Sumstats_for_GSEM_GWAS/Stroke_Gigastroke_2022_sumstats.gsem",
               "Data/Sumstats_for_GSEM_GWAS/PAD_2021_V2_sumstats.gsem",
               "Data/Sumstats_for_GSEM_GWAS/HF_2019_V2_sumstats.gsem")

ref= "Data/REF/reference.1000G.maf.0.005.txt"
CVD_names=c("CAD_2022",
            "Stroke_Gigastroke_2022",
            "PAD_2021_V2",
            "HF_2019_V2")

se.logit <- c(TRUE, TRUE, TRUE, TRUE)

CVD_sumstats <- sumstats(files = CVD_files,
                         ref = ref,
                         trait.names = CVD_names,
                         se.logit = se.logit)



saveRDS(CVD_sumstats, "Data/Sumstats_for_GSEM_GWAS/EFA_CVD_V2.5_sumstats.rds")