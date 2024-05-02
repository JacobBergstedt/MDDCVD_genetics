library(tidyverse)
library(glue)
library(stringr)
library(vroom)
library(fs)

group_keys <- c("Amygdala", "Anterior_cingulate_cortex", "Caudate", "Cerebellar_hemisphere", "Cerebellum", "Cortex", "Frontal_cortex", "Hippocampus", "Hypothalamus", "Nucleus_accumbens", "Putamen", "Spinal_cord", "Substantia_nigra")
loc <- map(glue("Data/REF/1000G_Phase3_baselineLD_v2.2_ldscores/baselineLD.{1:22}.annot.gz"), vroom) %>% 
  map(~ .[, 1:4])

concat_and_write_chr <- function(chr) {
  
  read_files <- function(i, chr) {
    
    annot <- read_delim(glue("Data/REF/GTEx_brain_1000Gv3_ldscores/GTEx_brain_1000Gv3_ldscores/GTEx_brain.{i}.{chr}.annot.gz"), col_types = "i", delim = "\t")
    ld <- read_delim(glue("Data/REF/GTEx_brain_1000Gv3_ldscores/GTEx_brain_1000Gv3_ldscores/GTEx_brain.{i}.{chr}.l2.ldscore.gz"), col_types = "ccid")
    M <- scan(glue("Data/REF/GTEx_brain_1000Gv3_ldscores/GTEx_brain_1000Gv3_ldscores/GTEx_brain.{i}.{chr}.l2.M"))
    M_5_50 <- scan(glue("Data/REF/GTEx_brain_1000Gv3_ldscores/GTEx_brain_1000Gv3_ldscores/GTEx_brain.{i}.{chr}.l2.M_5_50"))
    
    list(annot = annot,
         ld = ld,
         M = M,
         M_5_50 = M_5_50)  
    
  }
  
  flist <- map(1:13, read_files, chr = chr)
  annot_pref <- loc[[chr]]
  annot <- bind_cols(annot_pref, flist[[1]]$annot)
  ld <- flist[[1]]$ld
  M <- flist[[1]]$M
  M_5_50 <- flist[[1]]$M_5_50
  
  annot <- annot %>% dplyr::rename(Amygdala = ANNOT)
  ld <- ld %>% dplyr::rename(AmygdalaL2 = L2)
  
  for (i in 2:13) {
    
    annot_new <- flist[[i]]$annot %>% dplyr::rename("{group_keys[i]}" := ANNOT)
    ld_new <- flist[[i]]$ld %>% dplyr::rename("{group_keys[i]}L2" := L2)
    M_new <- flist[[i]]$M
    M_5_50_new <- flist[[i]]$M_5_50
    
    annot <- bind_cols(annot, annot_new)
    ld <- inner_join(ld, ld_new)
    M <- c(M, M_new)
    M_5_50 <- c(M_5_50, M_5_50_new)
    
  }
  
  write_delim(annot, glue("Data/REF/1000G_Phase3_GTEx_Brain/GTEx_brainLD.{chr}.annot.gz"), delim = "\t", quote = "none")
  write_delim(ld, glue("Data/REF/1000G_Phase3_GTEx_Brain/GTEx_brainLD.{chr}.l2.ldscore.gz"), delim = "\t", quote = "none")
  write(M, glue("Data/REF/1000G_Phase3_GTEx_Brain/GTEx_brainLD.{chr}.l2.M"), ncolumns = 13)
  write(M_5_50, glue("Data/REF/1000G_Phase3_GTEx_Brain/GTEx_brainLD.{chr}.l2.M_5_50"), ncolumns = 13)
  
}


map(1:22, concat_and_write_chr)