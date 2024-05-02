library(tidyverse)
library(LAVA)
loci <- read.loci("Data/REF/LAVA/blocks_s2500_m25_f1_w200.GRCh37_hg19.locfile")
loci_HLA <- loci %>% 
  filter(CHR == 6, START  > 28477797, STOP < 33448354) %>% 
  mutate(LOC = 1:16) %>% 
  add_row(LOC = 17, CHR = 6, START = 28477797, STOP = 33448354)

write_tsv(loci_HLA, "Data/REF/LAVA/HLA_loci.locfile")
