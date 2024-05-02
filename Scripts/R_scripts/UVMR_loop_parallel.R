library(TwoSampleMR)
library(dplyr)
source("./Scripts/R_scripts/Lib/functions_for_UVMR.R")
source("Scripts/R_scripts/Lib/functions_for_logistics.R")

CVD_factor_set <- get_MDD_var_list_updated_CVD()
MDD_set <- c("MDD_PGC3_V2","MDD_PGC3_noUKB_v2")





for (i in MDD_set) {
  for (j in CVD_factor_set) {
    UVMR(exposure = i,outcome = j, md_version = i)
  }
}
