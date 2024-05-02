genemat <- read_tsv("Data/Loci_gene_annotation.tsv")
snp_bim <- read.delim("/cluster/projects/p33/github/comorment/containers/reference/magma/g1000_eur/g1000_eur.bim", col.names = c("CHR", "RSID", "POSCM", "POS", "A1", "A2"), header = FALSE) %>% 
  as_tibble()
snp_bim_ranges <- GRanges(seqnames = paste0("chr", snp_bim$CHR), ranges = IRanges(start = snp_bim$POS, end = snp_bim$POS), RSID = snp_bim$RSID)


tissue <- genemat %>% select(CHR:STOP, Adipose_Tissue:Vagina)
tissue_names <- tissue %>% select(Adipose_Tissue:Vagina) %>% names()
tissue <- tissue %>% 
  rowwise() %>% 
  mutate(Top_tissue = tissue_names[which.max(c_across(Adipose_Tissue:Vagina))], Top_expr = max(c_across(Adipose_Tissue:Vagina))) %>% 
  ungroup()

tissue_gene_ranges <- tissue %>% 
  select(CHR:STOP, Top_tissue) %>%
  split(f = factor(.$Top_tissue, unique(.$Top_tissue))) %>% 
  map(~ GRanges(seqnames = .$CHR, ranges = IRanges(start = .$START,  end = .$STOP))) %>% 
  as("GRangesList")

overlaps <- findOverlaps(tissue_gene_ranges, snp_bim_ranges) %>% as_tibble() 
loci_top_tissue <- overlaps %>% 
  split(f = factor(.$queryHits, unique(.$queryHits))) %>% 
  map(~ snp_bim_ranges$RSID[.$subjectHits]) %>% 
  map_dfr(~ tibble(SNPS = paste0(., collapse = ";")), .id = "LOC")

write_tsv(loci_top_tissue, "Data/LAVA_loci_top_tissue_expression.tsv")
