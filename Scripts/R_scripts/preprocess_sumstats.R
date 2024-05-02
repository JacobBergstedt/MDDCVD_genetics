#!/cluster/software/EASYBUILD/R/4.1.0-foss-2021a/bin/Rscript
#SBATCH --job-name=preprocess_sumstats
#SBATCH --account=p33
#SBATCH --time=30:00
#SBATCH --cpus-per-task=12
#SBATCH --mem-per-cpu=10G

library(tidyverse)
library(fs)
library(vroom)
library(parallel)
source("./Scripts/R_scripts/Lib/functions_for_processing.R")
keys <- read_csv("Data/sumstats_keys_c.csv")
process_sumstats(c("Smoking"), keys, n_cores = 1)





