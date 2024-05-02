#!/cluster/software/EASYBUILD/R/4.1.0-foss-2021a/bin/Rscript
#SBATCH --job-name=LAVA_MDD
#SBATCH --output=./Slurm/Output/LAVA_local_cor_MDD_%A_%a.out
#SBATCH --error=./Slurm/Error/LAVA_local_cor_MDD_%A_%a.out
#SBATCH --account=p33
#SBATCH --time=7:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=25G
#SBATCH --array=1-100

library(LAVA)
library(tidyverse)
source("./Scripts/R_scripts/Lib/helpers.R")
source("./Scripts/R_scripts/Lib/functions_for_logistics.R")
source("./Scripts/R_scripts/Lib/functions_for_LAVA.R")

run_partial_correlation <- function(spec, locus, loc_info) {
  
  if (all(c(spec$target, spec$phenos) %in% locus$phenos)) {
    
    out <- try(run.pcor(locus, target = spec$target, phenos = spec$phenos, param.lim = 2, max.r2 = 0.99), silent = FALSE) 
    
    if (!grepl("error", attr(out, "class"))) out <- cbind(loc_info, out)
    
    out
    
  } else NULL
  
}
  

analyze_locus <- function(locus_pos, input) {
  
  locus <-  process.locus(locus = locus_pos, input = input)
  
  if ("MDD_Als_2023" %in% locus$phenos) {
    
    loc_info <-  data.frame(locus = locus$id, 
                            chr = locus$chr, 
                            start = locus$start, 
                            stop = locus$stop, 
                            n.snps = locus$n.snps, 
                            n.pcs = locus$K)  
    
    res_uni_locus <- run.univ(locus)
    res_uni <- bind_cols(loc_info, res_uni_locus)
    res_biv <- bind_cols(loc_info, run.bivar(locus, target = "MDD_Als_2023"))
      
    list(uni = res_uni, biv = res_biv)
    
  } else NULL
  
}



keys <- get_keys() %>%
  select(phenotype = Trait,
         cases = N_cases,
         controls =  N_controls,
         filename = Output_path) %>%
  filter(phenotype %in% c("MDD_Als_2023", "IL6"))

write.table(keys, file = "./Data/TMP_data/input_info_file_test_2024.txt", quote = FALSE)

input <- process.input(input.info.file = "./Data/TMP_data/input_info_file_test_2024.txt",
                       ref.prefix = "/cluster/projects/p33/github/comorment/magma/reference/magma/g1000_eur/g1000_eur",
                       phenos = keys$phenotype,
                       sample.overlap.file = "./Data/LDSC_sample_overlap_matrix_Als_2023.txt")

loci <- read.loci(paste0("./Data/TMP_data/LAVA_split/LAVA_loci_", 50, ".locfile")) %>%
  rowwise() %>% 
  dplyr::group_split()

locus <- loci[[1]]
res <- analyze_locus(locus, input)



