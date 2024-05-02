

# -------------------------------------------------------------------------



library(tidyverse)
library(ggplot2)
library(ggforestplot)
library(rlang)
source("Scripts/R_scripts/Lib/functions_for_logistics.R")

# 297 x 210 mm
# Margin = 25.4 mm

width <- 33.867
height <- 19.05
colors <- scales::brewer_pal(type = "qual", palette = 1)(8)[c(1, 2, 7, 4:6, 8)]
col_tib <- tibble(col = colors, 
                  Type = c("Cardiovascular disease", 
                           "Blood pressure", 
                           "Psychosocial", 
                           "Childhood maltreatment", 
                           "Metabolic", 
                           "Inflammation", 
                           "MDD")) %>% 
  mutate(Type = factor(Type, Type))
saveRDS(col_tib, "Data/color_keys.rds")


# Genome-wide correlation -------------------------------------------------
ldsc_correlations <- read_tsv("Results/LDSC_corr/ldsc_correlations_V2.tsv")

plt_frame <- ldsc_correlations %>% 
  select(-p1, -p2) %>% 
  filter(Trait1 == "MDD_PGC3_V2" | Trait2 == "MDD_PGC3_V2") %>% 
  mutate(Trait1 = ifelse(Trait1 != "MDD_PGC3_V2", Trait1, Trait2), Trait2 = "MDD_PGC3_V2") %>%
  filter(Trait1 %in% get_MDD_var_list()) %>%
  mutate(Type = case_when(
    Trait1 %in% c("Loneliness", "Overall_physical_activity_Doherty_2018_V2", "Sleep_duration_Doherty_2018_V2", "Educational_attainment", "High_physical_activity_Wang_2022_V2", "Smoking") ~ "Psychosocial",
    Trait1 %in% c("HDL_GLGC_2021", "LDL_GLGC_2021", "TC_GLGC_2021", "TG_GLGC_2021", "NonHDL_GLGC_2021", "BMI", "TIID_DIAMANTE_2022") ~ "Metabolic",
    Trait1 %in% c("CRP_CHARGEUKB_2022", "IL6") ~ "Inflammation",
    Trait1 %in% c("AF", "HF", "Stroke", "CAD", "PAD_2021_V2")  ~ "Cardiovascular disease",
    Trait1 %in% c("PP", "DBP", "SBP") ~ "Blood pressure",
    Trait1 %in% "CM_2021_V2" ~ "Childhood maltreatment"
  )
  ) %>% 
  mutate(Type = factor(Type, col_tib$Type)) %>% 
  group_by(Type) %>% 
  arrange(desc(abs(rg)), .by_group = TRUE) %>% 
  ungroup()


keys <- plt_frame %>%
  select(Trait = Trait1, Type) %>%
  bind_rows(tibble(Trait = "MDD_PGC3_V2", Type = "MDD")) %>% 
  mutate(Trait = factor(Trait, Trait), Type = factor(Type, levels(col_tib$Type))) %>% 
  mutate(Trait_labels = fct_recode(Trait,
                                   `Type II diabetes` = "TIID_DIAMANTE_2022",
                                   `Physical activity` = "Overall_physical_activity_Doherty_2018_V2",
                                   `Sleep duration` = "Sleep_duration_Doherty_2018_V2",
                                   `Educational attainment` = "Educational_attainment",
                                   `High physical activity` = "High_physical_activity_Wang_2022_V2",
                                   `Childhood maltreatment` = "CM_2021_V2",
                                   `Triglycerides` = "TG_GLGC_2021",
                                   `Total cholesterol` = "TC_GLGC_2021",
                                   `High-density lipo.` = "HDL_GLGC_2021",
                                   `Low-density lipo.` = "LDL_GLGC_2021",
                                   `Non high-density lipo.` = "NonHDL_GLGC_2021",
                                   `Peripheral artery dis.` = "PAD_2021_V2",
                                   `Heart failure` = "HF",
                                   `Coronary artery dis.` = "CAD",
                                   `Atrial fibrillation`= "AF",
                                   `Pulse pressure` = "PP",
                                   `Diastolic blood press.` = "DBP",
                                   `Systolic blood press.` = "SBP",
                                   `C-reactive protein` = "CRP_CHARGEUKB_2022",
                                   `MDD` = "MDD_PGC3_V2"))

saveRDS(keys, "Data/keys_types.rds")
saveRDS(col_tib, "Data/keys_color.rds")