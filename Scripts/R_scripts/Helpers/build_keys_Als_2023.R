
library(tidyverse)
library(ggplot2)
library(ggforestplot)
library(rlang)
source("Scripts/R_scripts/Lib/functions_for_logistics.R")

# 297 x 210 mm
# Margin = 25.4 mm

width <- 33.867
height <- 19.05
colors <- c("#7FC97F", "#BEAED4", "#BF5B17", "#E6DE47", "#386CB0", "#DE66A5", "#2D662D", "#BEAED4")
col_tib <- tibble(col = colors, 
                  Type = c("Cardiovascular dis.", 
                           "Blood pressure", 
                           "Psychosocial", 
                           "Childhood mal.", 
                           "Metabolic", 
                           "Inflammation", 
                           "MDD",
                           "MDD-ASCVD")) %>% 
  mutate(Type = factor(Type, Type))

# Genome-wide correlation -------------------------------------------------
ldsc_correlations <- read_tsv("Results/LDSC_corr/ldsc_correlations_Als_2023.tsv")

plt_frame <- ldsc_correlations %>% 
  select(-p1, -p2) %>% 
  filter(Trait1 == "MDD_Als_2023" | Trait2 == "MDD_Als_2023") %>% 
  mutate(Trait1 = ifelse(Trait1 != "MDD_Als_2023", Trait1, Trait2), Trait2 = "MDD_Als_2023") %>%
  filter(Trait1 %in% get_MDD_var_list_updated_CVD()) %>%
  mutate(Type = case_when(
    Trait1 %in% c("Loneliness", "Overall_physical_activity_Doherty_2018_V2", "Sleep_duration_Doherty_2018_V2", "Educational_attainment", "High_physical_activity_Wang_2022_V2", "Smoking") ~ "Psychosocial",
    Trait1 %in% c("HDL_GLGC_2021", "LDL_GLGC_2021", "TC_GLGC_2021", "TG_GLGC_2021", "NonHDL_GLGC_2021", "BMI", "TIID_DIAMANTE_2022") ~ "Metabolic",
    Trait1 %in% c("CRP_CHARGEUKB_2022", "IL6") ~ "Inflammation",
    Trait1 %in% c("AF", "HF_2019_V2", "Stroke_Gigastroke_2022", "CAD_2022", "PAD_2021_V2")  ~ "Cardiovascular dis.",
    Trait1 %in% c("PP", "DBP", "SBP") ~ "Blood pressure",
    Trait1 %in% "CM_2021_V2" ~ "Childhood mal."
  )
  ) %>% 
  mutate(Type = factor(Type, col_tib$Type)) %>% 
  group_by(Type) %>% 
  arrange(desc(abs(rg)), .by_group = TRUE) %>% 
  ungroup()


keys <- plt_frame %>%
  select(Trait = Trait1, Type) %>%
  bind_rows(tibble(Trait = "MDD_Als_2023", Type = "MDD")) %>% 
  mutate(Trait = factor(Trait, Trait), Type = factor(Type, levels(col_tib$Type))) %>% 
  mutate(Trait_labels = fct_recode(Trait,
                                   `Stroke` = "Stroke_Gigastroke_2022",
                                   `Type II diabetes` = "TIID_DIAMANTE_2022",
                                   `Physical activity` = "Overall_physical_activity_Doherty_2018_V2",
                                   `Sleep duration` = "Sleep_duration_Doherty_2018_V2",
                                   `Educational attainment` = "Educational_attainment",
                                   `Exercise` = "High_physical_activity_Wang_2022_V2",
                                   `Childhood maltreatment` = "CM_2021_V2",
                                   `Triglycerides` = "TG_GLGC_2021",
                                   `Total cholesterol` = "TC_GLGC_2021",
                                   `High-density lipo.` = "HDL_GLGC_2021",
                                   `Low-density lipo.` = "LDL_GLGC_2021",
                                   `Non high-density lipo.` = "NonHDL_GLGC_2021",
                                   `Peripheral artery dis.` = "PAD_2021_V2",
                                   `Heart failure` = "HF_2019_V2",
                                   `Coronary artery dis.` = "CAD_2022",
                                   `Atrial fibrillation`= "AF",
                                   `Pulse pressure` = "PP",
                                   `Diastolic blood press.` = "DBP",
                                   `Systolic blood press.` = "SBP",
                                   `C-reactive protein` = "CRP_CHARGEUKB_2022",
                                   `MDD` = "MDD_Als_2023")) %>% 
  mutate(Trait_label_short = fct_recode(Trait,
                                        `Stroke` = "Stroke_Gigastroke_2022",
                                        `T2D` = "TIID_DIAMANTE_2022",
                                        `Phys. act.` = "Overall_physical_activity_Doherty_2018_V2",
                                        `Sleep` = "Sleep_duration_Doherty_2018_V2",
                                        `Edu` = "Educational_attainment",
                                        `Exercise` = "High_physical_activity_Wang_2022_V2",
                                        `Child. mal.` = "CM_2021_V2",
                                        `TG` = "TG_GLGC_2021",
                                        `TC` = "TC_GLGC_2021",
                                        `HDL` = "HDL_GLGC_2021",
                                        `LDL` = "LDL_GLGC_2021",
                                        `NonHDL` = "NonHDL_GLGC_2021",
                                        `PAD` = "PAD_2021_V2",
                                        `HF` = "HF_2019_V2",
                                        `CAD` = "CAD_2022",
                                        `AF`= "AF",
                                        `PP` = "PP",
                                        `DBP` = "DBP",
                                        `SBP` = "SBP",
                                        `CRP` = "CRP_CHARGEUKB_2022",
                                        `MDD` = "MDD_Als_2023"))

saveRDS(keys, "Data/keys_types_Als_2023.rds")
saveRDS(col_tib, "Data/keys_color_Als_2023.rds")

