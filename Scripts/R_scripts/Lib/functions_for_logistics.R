get_keys <- function() {
  
  read_csv("./Data/sumstats_keys_c.csv")
  
}



get_ldsc_correlations <- function() {
  
  read_tsv("Results/LDSC_corr/ldsc_correlations_V2.tsv")
  
}


# get_MDD_var_list <- function() {
#   
#   c("AF",  
#     "CAD", 
#     "CRP_CHARGEUKB_2022", 
#     "Educational_attainment",
#     "HF_2019_V2_sumstats", 
#     "LDL_GLGC_2021", 
#     "PP", 
#     "Smoking", 
#     "DBP", 
#     "SBP", 
#     "HDL_GLGC_2021",
#     "TG_GLGC_2021",
#     "IL6", 
#     "PAD_2021_V2", 
#     "Stroke", 
#     "TC_GLGC_2021", 
#     "TIID_DIAMANTE_2022", 
#     "CM_2021_V2", 
#     "BMI",
#     "Sleep_duration_Doherty_2018_V2",
#     "Overall_physical_activity_Doherty_2018_V2",
#     "High_physical_activity_Wang_2022_V2",
#     "Loneliness",
#     "NonHDL_GLGC_2021")
#   
# }

get_MDD_var_list_without_CVD <- function() {
  
  c("CRP_CHARGEUKB_2022", 
    "Educational_attainment",
    "LDL_GLGC_2021", 
    "PP", 
    "Smoking",
    "DBP", 
    "SBP", 
    "HDL_GLGC_2021",
    "TG_GLGC_2021",
    "IL6", 
    "TC_GLGC_2021", 
    "TIID_DIAMANTE_2022", 
    "CM_2021_V2", 
    "BMI",
    "Sleep_duration_Doherty_2018_V2",
    "Overall_physical_activity_Doherty_2018_V2",
    "High_physical_activity_Wang_2022_V2",
    "Loneliness",
    "NonHDL_GLGC_2021")
  
}

get_MDD_var_list_only_updated_CVD <- function() {
  
  c("AF",  
    "CAD_2022",
    "HF_2019_V2",
    "PAD_2021_V2", 
    "Stroke_Gigastroke_2022")
  
  
}


get_MDD_var_list_updated_CVD <- function() {
  
  c(get_MDD_var_list_only_updated_CVD(),
    get_MDD_var_list_without_CVD())
  
}

get_MDD_var <- function() {
  
  "MDD_PGC3_V2"
  
}

get_MDD_Als <- function() {
  
  "MDD_Als_2023"
  
}

get_MDD_symptoms_list <- function() {
  
  c("Dep_symptoms_inadequacy",
    "Dep_symptoms_concentration",
    "Dep_symptoms_depressed",
    "Dep_symptoms_appetite",
    "Dep_symptoms_suicidal",
    "Dep_symptoms_anhedonia",
    "Dep_symptoms_sleep",
    "Dep_symptoms_psychomotor",
    "Dep_symptoms_tired")
  
}


get_ment_list <- function() {
  
  c("MDD_Als_2023_noUKB",
    "ADHD_IPsych_2023",
    "Insomnia_2019",
    "PTSD_2019",
    "SCZ_2022",
    "BIP_2021",
    "AlcDep_2018",
    "ANX_2020",
    "Neuroticism_2021",
    "CM_2021_V2")
  
}

