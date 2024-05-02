library(tidyverse)
library(glue)
source("Scripts/R_scripts/Lib/functions_for_logistics.R")
source("Scripts/R_scripts/Lib/helpers.R")
o <- get_keys() |> filter(Trait %in% c("MDD_Als_2023", get_MDD_var_list_updated_CVD()))
db <- glue("Data/Sumstats/{o$Trait}", "_sumstats") |> map(read_tsv)
names(db) <- o$Trait
sample_sizes <- db |> map_dfr(~ tibble(N = max(.$N), N_CAS = max(.$N_CAS), Prev = N_CAS / N), .id = "Trait") |> print(n = Inf)

paths <- dir_ls("Data/Sumstats")
paths_symptoms <- paths[grepl("[0-9]{5}", paths)]
symptoms <- paths_symptoms |> map(read_tsv)
names(symptoms) <- path_file(paths_symptoms)
sample_sizes <- symptoms |> map_dfr(~ tibble(N = max(.$N), N_CAS = max(.$N_CAS), Prev = N_CAS / N), .id = "Trait") |> print(n = Inf)
