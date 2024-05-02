#!/cluster/software/EASYBUILD/R/4.1.0-foss-2021a/bin/Rscript
#SBATCH --job-name=lhc_mr
#SBATCH --output=./Slurm/Output/lhc_mr_T2D.out
#SBATCH --error=./Slurm/Error/lhc_mr_T2D.out
#SBATCH --account=p33
#SBATCH --time=10:00:00
#SBATCH --cpus-per-task=6
#SBATCH --mem-per-cpu=25G

library(lhcMR)
library(tidyverse)
source("./Scripts/R_scripts/Lib/functions_for_logistics.R")
source("./Scripts/R_scripts/Lib/functions_for_lhc_mr.R")
select <- dplyr::select
slice <- dplyr::slice
#idx <- Sys.getenv("SLURM_ARRAY_TASK_ID")
idx <- 17
MDD <- read_tsv("Data/Sumstats/MDD_sumstats") %>% 
  filter(!CHR %in% c("Y", "X")) %>% 
  mutate(CHR = as.integer(CHR))


## File paths needed for the analysis
LD.filepath <-  "Data/REF/lhc_mr_required_ref/LDscores_filtered.csv" # LD scores
rho.filepath <-  "Data/REF/lhc_mr_required_ref/LD_GM2_2prm.csv" # local/SNP-specfic LD scores
ld <-  "Data/REF/lhc_mr_required_ref/eur_w_ld_chr/"
hm3 <-  "Data/REF/lhc_mr_required_ref/w_hm3.snplist"

var_list <- paste0(get_MDD_var_list(), "_sumstats")

###error: "BIP_sumstats","SCZ_sumstats"

outcome <- var_list[idx]
outcome_path <- paste("Data/Sumstats/", outcome, sep = "")
df_outcome <- read_tsv(outcome_path) %>% 
  filter(!CHR %in% c("Y", "X")) %>% 
  mutate(CHR = as.integer(CHR))

if (outcome %in% c("BMI_sumstats", "T2D_sumstats")) {
  
  df_outcome <- df_outcome %>% select(-Z)
  
}

trait.names=c("MDD", "df_outcome")
input.files = list(MDD, df_outcome)

df <- merge_sumstats(input.files, trait.names, LD.filepath, rho.filepath)

SP_list <-  calculate_SP(df, trait.names, run_ldsc = FALSE, run_MR = TRUE, hm3 = hm3, ld = ld, nStep = 2, SP_single = 3, SP_pair = 50, SNP_filter = 10, nCores = 5)
res <- lhc_mr(SP_list, trait.names, paral_method = "lapply", nCores = 5, nBlock = 200)
res <- as.data.frame(res)
write_tsv(res, paste0("Results/Chunk_data/Lhc_MR/", outcome, ".tsv"))

