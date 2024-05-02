library(GenomicRanges)
gene_mat <- read_tsv("Data/geneMatrix.tsv.gz") %>% 
  filter(gene_type == "protein_coding", !is.na(hg19g0))

no_tissue_expression <- gene_mat %>% 
  select(Adipose_Tissue:Vagina) %>% 
  rowwise() %>% 
  summarize(nr_na = sum(is.na(c_across(Adipose_Tissue:Vagina))))

gene_mat <- gene_mat %>% filter(no_tissue_expression$nr_na < 20)

gene_mat_sub <- gene_mat %>% filter(!grepl("chrM|chrY|chrX", hg19g0))
gene_ranges <- GRanges(seqnames = gene_mat_sub$hg19g0, 
                       IRanges(start = gene_mat_sub$g1, end = gene_mat_sub$g2), 
                       mcols = DataFrame(ID = gene_mat_sub$gene_id))
loci <- read.loci("/cluster/projects/p33/github/comorment/containers/reference/lava/blocks_s2500_m25_f1_w200.GRCh37_hg19.locfile")
whole_genome_tib <- as_tibble(loci) %>% group_by(CHR) %>% summarize(START = min(START), STOP = max(STOP))
whole_genome_ranges <- GRanges(seqnames = paste0("chr", whole_genome_tib$CHR), 
                               ranges = IRanges(start = whole_genome_tib$START, end = whole_genome_tib$STOP))

missing_from_gene_mat <- setdiff(whole_genome_ranges, gene_ranges)
missing_from_gene_mat_narrow <- missing_from_gene_mat[width(ranges(missing_from_gene_mat)) < 50]

for (i in seq_along(missing_from_gene_mat_narrow)) {
  
  nearest_gene_seq <- nearest(missing_from_gene_mat_narrow[i], gene_ranges)
  expanded_range <- union(gene_ranges[nearest_gene_seq], missing_from_gene_mat_narrow[i])
  mcols(expanded_range) <- mcols(gene_ranges[nearest_gene_seq])
  gene_ranges[nearest_gene_seq] <- expanded_range
  
}

missing_from_gene_mat <- setdiff(whole_genome_ranges, gene_ranges)
missing_from_gene_mat_tiled <- missing_from_gene_mat %>% tile(n = 2) %>% unlist()
for (i in seq_along(missing_from_gene_mat)) {
  
  if (i%%1000 == 0) print(i)
  nearest_gene_seq <- nearest(missing_from_gene_mat[i], gene_ranges)
  expanded_range <- union(gene_ranges[nearest_gene_seq], missing_from_gene_mat[i])
  mcols(expanded_range) <- mcols(gene_ranges[nearest_gene_seq])
  gene_ranges[nearest_gene_seq] <- expanded_range
  
}

gene_ranges_tib <- gene_ranges %>% 
  as_tibble() %>% 
  dplyr::rename(CHR = seqnames, START = start, STOP = end, gene_id = mcols.ID) %>% 
  inner_join(gene_mat_sub) %>% 
  select(-starts_with("hg"), -g2, -g1, -h1, -h2)

write_tsv(gene_ranges_tib, "Data/Loci_gene_annotation.tsv")

