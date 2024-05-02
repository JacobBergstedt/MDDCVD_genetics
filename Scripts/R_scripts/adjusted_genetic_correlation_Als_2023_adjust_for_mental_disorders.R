
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
  
  res <- m[grepl("T1~~T2", paste0(m$lhs, m$op, m$rhs)), c("STD_Genotype" , "STD_Genotype_SE", "p_value")]
  tibble(rg = res[[1]], rg_se = as.numeric(res[[2]]), P = as.numeric(res[[3]]))
  
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


fit_adjusted_model <- function(trait1, trait2, adjust) {
  
  res_traits <- tibble(Trait1 = trait1, Trait2 = trait2)
  spec <- glue("T1 =~ {trait1}", 
               "T2 =~ {trait2}", 
               "T1 ~ {adjust}", 
               "T2 ~ {adjust}", 
               .sep = "\n")
  
  m <- usermodel(gen_cov, 
                 estimation = "DWLS",
                 model = spec, 
                 CFIcalc = TRUE, 
                 std.lv = TRUE,
                 imp_cov = FALSE)
  
  res <- get_cor(m$results)
  res$Adjust_for <- adjust
  bind_cols(res_traits, res)
  
  
  
}

fit_models <- function(trait1, trait2, adjustment) {
  
  if (adjustment == "Nothing") res <- fit_unadjusted_model(trait1, trait2)
  else res <- fit_adjusted_model(trait1, trait2, adjust = adjustment)
  res
  
}


# mental_disorder_traits <- c("MDD", "MDD_Howard", "ANX", "SCZ", "BIP", "ANX + MDD_Howard")
cvd_traits <- c("CAD_2022 + HF_2019_V2 + Stroke_Gigastroke_2022 + PAD_2021_V2",
                get_MDD_var_list_only_updated_CVD())

# Main --------------------------------------------------------------------
# keys_types <- readRDS("Data/keys_types_Als_2023.rds") %>% filter(Trait %in% get_MDD_var_list_without_CVD())
# keys_label <- keys_types %>% select(Trait, Trait_labels)
# 
# type_models <- keys_types %>% 
#   filter(!Trait %in% c("PP", "SBP", "DBP", "CM_2021_V2")) %>% 
#   group_by(Type) %>% 
#   summarize(Model = paste0(Trait, collapse = " + ")) %>% 
#   ungroup() %>% 
#   pull(Model)

gen_cov <- readRDS("Data/Genetic_covariance_matrices/genetic_covariance_Als_2023_mental_disorders.rds")

grid <- expand.grid(trait1 = "MDD_Als_2023", 
                    trait2 = cvd_traits, 
                    adjustment = c("Nothing", "PTSD_2019", "BIP_2021", "SCZ_2022", "ANX_2020"), stringsAsFactors = FALSE)

res <- pmap_dfr(grid, fit_models)

keys_label <- bind_rows(keys_label,
                        tibble(Trait = c("Nothing", 
                                         "Loneliness + Smoking + High_physical_activity_Wang_2022_V2 + Educational_attainment + Overall_physical_activity_Doherty_2018_V2 + Sleep_duration_Doherty_2018_V2",
                                         "TIID_DIAMANTE_2022 + TG_GLGC_2021 + HDL_GLGC_2021 + BMI + NonHDL_GLGC_2021 + TC_GLGC_2021 + LDL_GLGC_2021",
                                         "IL6 + CRP_CHARGEUKB_2022"),
                               Trait_labels = c("Nothing", "Psychosocial", "Metabolic", "Inflammation"))) %>% 
  rename(Adjust_for = Trait, Adjust_for_label = Trait_labels)

res_out <- res %>%
  mutate(Adjust_for = fct_rev(factor(Adjust_for, unique(Adjust_for)))) %>%
  mutate(Trait2 = recode(Trait2, `CAD + HF + Stroke + PAD_2021_V2` = "CVD")) %>%
  inner_join(keys_label)

saveRDS(res_out, "Results/MDD_CVD_adjusted_genetic_correlations_Als_2023.rds")
saveRDS(res_out, paste0(Sys.getenv("export"), "/Adjusted_gencor_Als_2023.rds"))



