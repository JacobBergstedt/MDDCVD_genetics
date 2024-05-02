library(jsonlite)
library(tidyverse)
library(glue)
source("Scripts/R_scripts/Lib/functions_for_logistics.R")
traits <- c("MDD_Als_2023", get_MDD_var_list_updated_CVD())
paths <- map(traits, ~ glue("Results/Mixer/MDD_Als_2023/{.}.fit.rep{1:20}.json"))

res <- vector(mode = "list", length = length(traits))
names(res) <- traits

for(i in seq_along(paths)) {
  
  p <- map(paths[[i]], ~ read_json(.)$ci)
  res[[i]] <- map_dfr(p, unlist, .id = "Rep")
  
}

res <- map_dfr(res, ~ ., .id = "Trait")

write_tsv(res, "Results/res_Mixer_uni_repeats.tsv")


# -------------------------------------------------------------------------

traits <- c("MDD_Als_2023", 
            "PAD",
            get_MDD_var_list_updated_CVD())

paths_fit <- map(traits, ~ glue("Results/Mixer/MDD_Als_2023/Output/{.}.fit.output.csv"))
paths_test <- map(traits, ~ glue("Results/Mixer/MDD_Als_2023/Output/{.}.test.output.csv"))

res_uni_fit <- map_dfr(paths_fit, read_delim)
res_uni_test <- map_dfr(paths_test, read_delim)

write_tsv(res_uni_fit, "Results/res_mixer_univariate_fit.tsv")
write_tsv(res_uni_test, "Results/res_mixer_univariate_test.tsv")