library(GenomicSEM)
library(tidyverse)
library(data.table)
library(glue)
library(Matrix)
library(stats)
library(lavaan)

# Trait information
source("Scripts/R_scripts/Lib/functions_for_logistics.R")
keys <- get_keys()
var_list <- get_MDD_var_list()
keys %>% filter(Trait %in% var_list)

##################################################################################################
###Compare mediation versus covariation models for CVD and MDD 

LDSCoutput <- readRDS("Data/Genetic_covariance_matrices/genetic_covariance_Als_2023.rds")

### 'Mediation' models

fmediation <- function(trait) {
  mediationmodel <- glue("CVD =~ 1*CAD_2022 + HF_2019_V2 + PAD_2021_V2 + Stroke_Gigastroke_2022",
                         "CVD ~ MDD_Als_2023 + {trait}",                     
                         "{trait} ~ MDD_Als_2023",
                         .sep = "\n")
  m <- usermodel(covstruc = LDSCoutput, estimation="DWLS", model=mediationmodel)
  med_res <- tibble(m$results)
  med_res$trait <- {trait}
  as_tibble(med_res)
}


traits <- c("Smoking", "PP", "DBP", "SBP", "CRP_CHARGEUKB_2022", "Educational_attainment", "TC_GLGC_2021", "TG_GLGC_2021","LDL_GLGC_2021", "HDL_GLGC_2021", "NonHDL_GLGC_2021", "CM_2021_V2", "IL6", 
            "BMI", "TIID_DIAMANTE_2022", "Overall_physical_activity_Doherty_2018_V2", "High_physical_activity_Wang_2022_V2","Sleep_duration_Doherty_2018_V2", "Loneliness",
            "DBP + SBP + PP", "TG_GLGC_2021 + TC_GLGC_2021 + LDL_GLGC_2021 + HDL_GLGC_2021 +  NonHDL_GLGC_2021 + BMI + TIID_DIAMANTE_2022", "IL6 + CRP_CHARGEUKB_2022", 
            "Smoking + Educational_attainment + Overall_physical_activity_Doherty_2018_V2 + High_physical_activity_Wang_2022_V2 + Sleep_duration_Doherty_2018_V2 + Loneliness + CM_2021_V2"
)

mediation <- traits %>% 
  map(fmediation)

names(mediation) <- traits

### Covariation models


fcovariation <- function(trait) {
  covarmodel <- glue("CVD =~ 1*CAD_2022 + HF_2019_V2 + PAD_2021_V2 + Stroke_Gigastroke_2022",
                     "CVD ~ MDD_Als_2023",                     
                     "{trait} ~ MDD_Als_2023",
                     "{trait} ~~ CVD",
                     .sep = "\n")
  m <- usermodel(covstruc = LDSCoutput, estimation="DWLS", model=covarmodel)
  cov_res <- tibble(m$results)
  cov_res$trait <- {trait}
  as_tibble(cov_res)
}

covariation <- traits %>% 
  map(fcovariation)

names(covariation) <- traits

### Compare

comp <- function(trait) {
  
  med <- mediation[[trait]] %>% filter(lhs=="CVD", rhs=="MDD_Als_2023")
  cov <- covariation[[trait]] %>% filter(lhs=="CVD", rhs=="MDD_Als_2023")
  lCI_med <- med$STD_Genotype - (1.96*as.numeric(med$STD_Genotype_SE))
  uCI_med <- med$STD_Genotype + (1.96*as.numeric(med$STD_Genotype_SE))
  lCI_cov <- cov$STD_Genotype - (1.96*as.numeric(cov$STD_Genotype_SE))
  uCI_cov <- cov$STD_Genotype + (1.96*as.numeric(cov$STD_Genotype_SE))
  attenuation <- (lCI_cov - lCI_med)/ lCI_cov *100
  sig <- if(lCI_cov>uCI_med | uCI_cov<lCI_med) {"yes"} else {"no"}
  tibble(lCI_med, uCI_med, lCI_cov, uCI_cov, attenuation, sig)
}

dif <- traits %>%
  map(comp)

names(dif) <- traits

res <- rbindlist(dif, idcol=T)
fwrite(res, "Results/gSEM_joelle/results_mediation_vs_covariation_Als.tsv")


