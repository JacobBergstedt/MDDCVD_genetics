#!/cluster/software/EASYBUILD/R/4.1.0-foss-2021a/bin/Rscript
#SBATCH --job-name=GSEM_GWAS_MDDCVD
#SBATCH --output=./Slurm/Output/GSEM_GWAS_MDD_CVD_%A_%a.out
#SBATCH --error=./Slurm/Error/GSEM_GWAS_MDD_CVD_%A_%a.out
#SBATCH --account=p33_tsd
#SBATCH --time=4:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=7G
#SBATCH --array=1-300

library(GenomicSEM)
library(tidyverse)
library(lavaan)
library(glue)

IDX <- Sys.getenv("SLURM_ARRAY_TASK_ID")
gen_cov <- readRDS("Data/Genetic_covariance_matrices/genetic_covariance_for_GWAS_MDDCVD_Als_2023.rds")
snps <- readRDS(glue("Data/TMP_data/GSEM_split/MDD_CVD_GSEM_LOCI_SPLIT{IDX}.rds"))

model_common <- '
CVD =~ 1*CAD_2022 + Stroke_Gigastroke_2022 + PAD_2021_V2 + HF_2019_V2
CVDMDD =~ 1*MDD_Als_2023 + CVD
CVD ~~ 0*CVD
CVDMDD ~ SNP
'

common_pathway <- userGWAS(covstruc = gen_cov, 
                           SNPs = snps, 
                           estimation = "DWLS",
                           model = model_common, 
                           printwarn = TRUE,
                           cores = 1, 
                           toler = FALSE,
                           SNPSE = FALSE, 
                           parallel = FALSE)

saveRDS(common_pathway, paste0("Results/Chunk_data/GSEM_MDDCVD/Als_2023/GSEM_MDD_Als_2023_CVD_V2.5_COMMON_PATHWAY_", IDX, ".rds"))

model_independent <- '
CVD =~ 1*CAD_2022 + Stroke_Gigastroke_2022 + PAD_2021_V2 + HF_2019_V2
CVDMDD =~ 1*MDD_Als_2023 + CVD
CVD ~~ 0*CVD
CAD_2022 + Stroke_Gigastroke_2022 + PAD_2021_V2 + HF_2019_V2 + MDD_Als_2023 ~ SNP
'

independent_pathway <- userGWAS(covstruc = gen_cov, 
                                SNPs = snps, 
                                estimation = "DWLS",
                                model = model_independent, 
                                printwarn = TRUE,
                                cores = 1, 
                                toler = FALSE,
                                SNPSE = FALSE, 
                                parallel = FALSE)


saveRDS(independent_pathway, paste0("Results/Chunk_data/GSEM_MDDCVD/Als_2023/GSEM_MDD_Als_2023_CVD_V2.5_INDEPENDENT_PATHWAY_", IDX, ".rds"))

