#!/cluster/software/EASYBUILD/R/4.1.0-foss-2021a/bin/Rscript
#SBATCH --job-name=gSEM_LDSC
#SBATCH --output=./Slurm/Output/gSEM_LDSC.out
#SBATCH --error=./Slurm/Error/gSEM_LDSC.out
#SBATCH --account=p33
#SBATCH --time=1:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=7G

library(tidyverse)
library(GenomicSEM)
library(glue)
source("./Scripts/R_scripts/Lib/functions_for_logistics.R")


# Main --------------------------------------------------------------------


REF <- "/cluster/projects/p33/github/comorment/containers/reference/ldsc/"
var_list <- c("MDD_Als_2023", get_MDD_var_list_updated_CVD())

keys <- get_keys() %>% 
  mutate(Output_munged = paste0(gsub("/Sumstats/", "/Munged_sumstats/", Output_path), "_munged.sumstats"),
         Sample_prev = N_cases / (N_cases + N_controls)) %>% 
  filter(Trait %in% var_list)

ld <- "Data/REF/eur_w_ld_chr"
wld <- ld
gen_cov <- ldsc(keys$Output_munged, 
                sample.prev = keys$Sample_prev, 
                population.prev = keys$Pop_prev, 
                ld = ld, 
                wld = wld, 
                trait.names = keys$Trait)

saveRDS(gen_cov, "Data/Genetic_covariance_matrices/genetic_covariance_Als_2023.rds")
