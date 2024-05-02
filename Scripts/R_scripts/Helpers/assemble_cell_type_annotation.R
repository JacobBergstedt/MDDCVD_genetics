library(tidyverse)
library(glue)
library(stringr)

group_keys <- read_tsv("Data/REF/1000G_Phase3_cell_type_groups/1000G_Phase3_cell_type_groups/names") %>% 
  deframe()

group_keys[] <- str_sub(group_keys, end = -5)


concat_and_write_chr <- function(chr) {
  
  read_files <- function(i, chr) {
    
    annot <- read_delim(glue("Data/REF/1000G_Phase3_cell_type_groups/1000G_Phase3_cell_type_groups/cell_type_group.{i}.{chr}.annot.gz"), col_types = "cicdd")
    ld <- read_delim(glue("Data/REF/1000G_Phase3_cell_type_groups/1000G_Phase3_cell_type_groups/cell_type_group.{i}.{chr}.l2.ldscore.gz"), col_types = "ccid")
    M <- scan(glue("Data/REF/1000G_Phase3_cell_type_groups/1000G_Phase3_cell_type_groups/cell_type_group.{i}.{chr}.l2.M"))
    M_5_50 <- scan(glue("Data/REF/1000G_Phase3_cell_type_groups/1000G_Phase3_cell_type_groups/cell_type_group.{i}.{chr}.l2.M_5_50"))
    
    list(annot = annot,
         ld = ld,
         M = M,
         M_5_50 = M_5_50)
    
  }
  
  
  flist <- map(1:10, read_files, chr = chr)
  
  annot <- flist[[1]]$annot
  ld <- flist[[1]]$ld
  M <- flist[[1]]$M
  M_5_50 <- flist[[1]]$M_5_50
  
  annot <- annot %>% dplyr::rename(Adrenal_pancreas = ANNOT)
  ld <- ld %>% dplyr::rename(Adrenal_pancreasL2 = L2)
  
  for (i in 2:10) {
    
    
    annot_new <- flist[[i]]$annot %>% dplyr::rename("{group_keys[i]}" := ANNOT) %>% select(-CM)
    ld_new <- flist[[i]]$ld %>% dplyr::rename("{group_keys[i]}L2" := L2)
    M_new <- flist[[i]]$M
    M_5_50_new <- flist[[i]]$M_5_50
    
    annot <- inner_join(annot, annot_new)
    ld <- inner_join(ld, ld_new)
    M <- c(M, M_new)
    M_5_50 <- c(M_5_50, M_5_50_new)
    
  }
  
  write_delim(annot, glue("Data/REF/1000G_Phase3_cell_type_groups_concat/cell_type_groupsLD.{chr}.annot.gz"), delim = "\t", quote = "none")
  write_delim(ld, glue("Data/REF/1000G_Phase3_cell_type_groups_concat/cell_type_groupsLD.{chr}.l2.ldscore.gz"), delim = "\t", quote = "none")
  write(M, glue("Data/REF/1000G_Phase3_cell_type_groups_concat/cell_type_groupsLD.{chr}.l2.M"), ncolumns = 10)
  write(M_5_50, glue("Data/REF/1000G_Phase3_cell_type_groups_concat/cell_type_groupsLD.{chr}.l2.M_5_50"), ncolumns = 10)
  
}


 map(1:22, concat_and_write_chr)