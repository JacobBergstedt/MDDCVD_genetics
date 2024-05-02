#!/cluster/software/EASYBUILD/R/4.1.0-foss-2021a/bin/Rscript
#SBATCH --job-name=LAVA_MDD_Als_2023
#SBATCH --output=./Slurm/Output/LAVA_local_cor_MDD_Als_2023_%A_%a.out
#SBATCH --error=./Slurm/Error/LAVA_local_cor_MDD_Als_2023_%A_%a.out
#SBATCH --account=p33
#SBATCH --time=4:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=25G
#SBATCH --array=1-3

library(LAVA)
library(tidyverse)
library("optparse")
source("./Scripts/R_scripts/Lib/helpers.R")
source("./Scripts/R_scripts/Lib/functions_for_logistics.R")
source("./Scripts/R_scripts/Lib/functions_for_LAVA.R")


option_list <- list( make_option("--trait", type="character"))
opt <-  parse_args(OptionParser(option_list = option_list))

target_trait <- "MDD_Als_2023"
CVD_trait <- opt$trait

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

idx_split <- Sys.getenv("SLURM_ARRAY_TASK_ID")
keys <- get_keys() %>%
  select(phenotype = Trait,
         cases = N_cases,
         controls =  N_controls,
         filename = Output_path) %>%
  filter(phenotype %in% c("MDD_Als_2023", CVD_trait))

keys <- keys %>%
  mutate(filename = if_else(phenotype == "MDD_Als_2023", "Data/TMP_data/MDD_Als_2023_sumstats_HM3_incl_HLA.tsv", filename))

if (idx_split == 1) {

  write.table(keys, file = "./Data/TMP_data/input_info_file.txt", quote = FALSE)

}

input <- process.input(input.info.file = "./Data/TMP_data/input_info_file.txt",
                       ref.prefix = "/cluster/projects/p33/github/comorment/magma/reference/magma/g1000_eur/g1000_eur",
                       phenos = keys$phenotype,
                       sample.overlap.file = "./Data/LDSC_sample_overlap_matrix_Als_2023.txt")

loci <- read.loci(paste0("./Data/TMP_data/LAVA_split/LAVA_loci_", idx_split, ".locfile")) %>%
  rowwise() %>%
  dplyr::group_split()

res <- map(loci, analyze_locus, input = input)
res_uni <- compact(res) %>% map("uni") %>% bind_rows()
res_biv <- compact(res) %>% map("biv") %>% bind_rows()

write_tsv(res_uni, paste0("./Results/Chunk_data/LAVA_HapMap3/LAVA_uni_", CVD_trait, "_",  idx_split, ".tsv"))
write_tsv(res_biv, paste0("./Results/Chunk_data/LAVA_HapMap3/LAVA_biv_", CVD_trait, "_", idx_split, ".tsv"))


