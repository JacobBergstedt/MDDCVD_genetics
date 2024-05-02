#!/cluster/software/EASYBUILD/R/4.1.0-foss-2021a/bin/Rscript
#SBATCH --job-name=gSEM_LDSC
# SBATCH --output=./Slurm/Output/gSEM_LDSC.out
#SBATCH --error=./Slurm/Error/gSEM_LDSC.out
#SBATCH --account=p33
#SBATCH --time=5:00:00
#SBATCH --cpus-per-task=2
#SBATCH --mem-per-cpu=30G

library(tidyverse)
library(GenomicSEM)
library(glue)
source("./Scripts/R_scripts/Lib/functions_for_logistics.R")

# Main --------------------------------------------------------------------

REF <- "/cluster/projects/p33/github/comorment/ldsc/reference/"
keys <- get_keys() %>% 
  mutate(Output_munged = paste0(gsub("/Sumstats/", "/Munged_sumstats/", Output_path), "_munged.sumstats"),
         Sample_prev = N_cases / (N_cases + N_controls))


var_list <- c(get_MDD_Als(), get_MDD_var_list_only_updated_CVD())[-2]
keys_MDD_CVD <- keys %>% 
  filter(Trait %in% var_list) %>% 
  select(Output_munged, Sample_prev, Pop_prev, Trait)

ld <- "Data/REF/eur_w_ld_chr"
wld <- ld
gen_cov <- ldsc(keys_MDD_CVD$Output_munged, 
                sample.prev = keys_MDD_CVD$Sample_prev, 
                population.prev = keys_MDD_CVD$Pop_prev, 
                ld = ld, 
                wld = wld, 
                ldsc.log = "ldsc.log",
                trait.names = keys_MDD_CVD$Trait)

saveRDS(gen_cov, "Data/Genetic_covariance_matrices/genetic_covariance_for_GWAS_MDDCVD_Als_2023.rds")


var_list <- get_MDD_var_list_only_updated_CVD()[-1]
keys_CVD <- keys %>% 
  filter(Trait %in% var_list) %>% 
  select(Output_munged, Sample_prev, Pop_prev, Trait)

ld <- "Data/REF/eur_w_ld_chr"
wld <- ld
gen_cov <- ldsc(keys_CVD$Output_munged, 
                sample.prev = keys_CVD$Sample_prev, 
                population.prev = keys_CVD$Pop_prev, 
                ld = ld, 
                wld = wld, 
                ldsc.log = "ldsc.log",
                trait.names = keys_CVD$Trait)

saveRDS(gen_cov, "Data/Genetic_covariance_matrices/genetic_covariance_for_GWAS_CVD_V2.5.rds")



