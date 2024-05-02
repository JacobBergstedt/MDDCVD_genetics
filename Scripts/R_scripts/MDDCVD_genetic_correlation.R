
#!/cluster/software/EASYBUILD/R/4.1.0-foss-2021a/bin/Rscript
#SBATCH --job-name=adjusted_genetic_correlation
#SBATCH --output=./Slurm/Output/adjusted_genetic_correlation.out
#SBATCH --error=./Slurm/Error/adjusted_genetic_correlation.out
#SBATCH --account=p33
#SBATCH --time=1:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=20G


library(tidyverse)
library(GenomicSEM)
library(glue)
source("./Scripts/R_scripts/Lib/functions_for_logistics.R")

get_cor <- function(m) {
  
  res <- m[grepl("T1~~T2", paste0(m$lhs, m$op, m$rhs)), c("STD_Genotype" , "STD_Genotype_SE")]
  tibble(rg = res[[1]], rg_se = as.numeric(res[[2]]))
  
}



fit_unadjusted_model <- function(trait1, trait2) {
  
  res_traits <- tibble(Trait1 = trait1, Trait2 = trait2)
  spec <- glue("T1 =~ {trait1}", "T2 =~ {trait2}", .sep = "\n")
  m <- usermodel(gen_cov, 
                 estimation = "DWLS",
                 model = spec, 
                 CFIcalc = TRUE, 
                 std.lv = TRUE,
                 imp_cov = FALSE)
  res <- get_cor(m$results)
  res$Adjust_for <- "Nothing"
  bind_cols(res_traits, res)
  
  
}


fit_MDDCVD_model <- function(trait) {
  
  spec <- glue("CVD =~ CAD + HF + PAD + Stroke", 
               "CVDMDD =~ 1 * CVD + MDD", 
               "CVD ~~ 0 * CVD",
               "CVDMDD ~~ {trait}",
               .sep = "\n")
  
  
  m <- usermodel(gen_cov, 
                 estimation = "DWLS",
                 model = spec, 
                 CFIcalc = TRUE, 
                 std.lv = TRUE,
                 imp_cov = FALSE)
  print(m$modelfit)
  res <- m$results
  as_tibble(res[grepl(glue("CVDMDD~~{trait}"), paste0(res$lhs, res$op, res$rhs)), c("STD_Genotype" , "STD_Genotype_SE")])
  
}
  
# Main --------------------------------------------------------------------

traits <- c("CRP", "Educational_attainment", "Smoking", "BMI", "LDL", "TC", "TG", "T2D", "PP", "DBP", "SBP", "HDL", "IL6")
gen_cov <- readRDS("Data/Genetic_covariance_matrices/genetic_covariance.rds")
res <- map_dfr(set_names(traits, traits), fit_MDDCVD_model, .id = "Trait")
saveRDS(res, "Results/MDDCVD_latent_factor_genetic_correlations.rds")



