#!/cluster/software/EASYBUILD/R/4.1.0-foss-2021a/bin/Rscript
#SBATCH --job-name=functional_enrichment
#SBATCH --output=./Slurm/Output/functional_enrichments.out
#SBATCH --error=./Slurm/Error/functional_enrichments.out
#SBATCH --account=p33
#SBATCH --time=2:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=50G

library(tidyverse)
library(GenomicSEM)
library(glue)
source("./Scripts/R_scripts/Lib/functions_for_logistics.R")
source("./Scripts/R_scripts/Lib/helpers.R")

# gen_cov_unstrat <- readRDS("Data/Genetic_covariance_matrices/genetic_covariance.rds")
gen_cov <- readRDS("Data/Genetic_covariance_matrices/stratified_genetic_covariance_V2.5.rds")
# gen_cov_brain <- readRDS("Data/Genetic_covariance_matrices/stratified_genetic_covariance_baseline_brain.rds")
# gen_cov_tissues <- readRDS("Data/Genetic_covariance_matrices/stratified_genetic_covariance_baseline_tissues.rds")
# tissues <- paste0(gen_cov_tissues$Select[97:106, 1], "L2")
# brain_tissues <- paste0(gen_cov_brain$Select$V1[97:109], "L2")


spec <- 'CVD =~ 1*CAD_2022 + HF_2019_V2 + PAD_2021_V2 + Stroke_Gigastroke_2022
CVDMDD =~ 1*MDD_PGC3_V2 + CVD
CVD ~~ 0*CVD
'

params <- "CVDMDD ~~ CVDMDD"
res <- as_tibble(enrich(s_covstruc = gen_cov, model = spec, std.lv = TRUE, params = params)[[1]]) %>% 
  arrange(desc(Enrichment - 1.96 * Enrichment_SE))

saveRDS(res, "Results/MDDCVD_enrichment_baseline_V2.5.rds")
saveRDS(res, "/tsd/p33/data/durable/file-export/MDDCVD_enrichment_baseline_V2.5.rds")