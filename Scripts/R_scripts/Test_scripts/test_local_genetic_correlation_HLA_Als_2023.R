library(LAVA)
library(tidyverse)
source("./Scripts/R_scripts/Lib/helpers.R")
source("./Scripts/R_scripts/Lib/functions_for_logistics.R")
source("./Scripts/R_scripts/Lib/functions_for_LAVA.R")

analyze_locus <- function(locus_pos, input) {
  
  locus <-  process.locus(locus = locus_pos, input = input)
  if (!is_null(locus)) {
    loc_info <-  data.frame(locus = locus$id, chr = locus$chr, start = locus$start, stop = locus$stop, n.snps = locus$n.snps, n.pcs = locus$K)  
    res_uni_locus <- run.univ(locus)
    res_uni <- bind_cols(loc_info, res_uni_locus)
    res_biv <- bind_cols(loc_info, run.bivar(locus, param.lim = 10))
    list(uni = res_uni, biv = res_biv)
  } else NULL
  
}

loci <- read.loci("./Data/REF/LAVA/HLA_loci.locfile") %>%
  rowwise() %>%
  dplyr::group_split()

res <- vector(mode = "list", length = 5)
res_CVD_traits <- get_MDD_var_list_only_updated_CVD()
names(res) <- res_CVD_traits

for (trait in res_CVD_traits) {
  
  keys <- get_keys() %>%
    select(phenotype = Trait,
           cases = N_cases,
           controls =  N_controls,
           filename = Output_path) %>%
    filter(phenotype %in% c(trait, "MDD_Als_2023"))
  
  keys <- keys %>%
    mutate(filename = if_else(phenotype == "MDD_Als_2023", "Data/TMP_data/MDD_Als_2023_sumstats_HM3_incl_HLA.tsv", filename))
  
  path_input_info_file <- paste0("./Data/TMP_data/input_info_file_", trait, ".txt")
  
  write.table(keys, file = path_input_info_file, quote = FALSE)
  
  input <- process.input(input.info.file = path_input_info_file,
                         ref.prefix = "/cluster/projects/p33/github/comorment/magma/reference/magma/g1000_eur/g1000_eur",
                         phenos = keys$phenotype,
                         sample.overlap.file = "./Data/LDSC_sample_overlap_matrix_Als_2023.txt")
  
  
  res[[trait]] <- map(loci, analyze_locus, input)
  
}






# "Data/TMP_data/MDD_Als_2023_sumstats_HM3_incl_HLA.tsv"