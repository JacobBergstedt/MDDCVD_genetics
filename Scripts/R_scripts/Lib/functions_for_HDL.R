

munge_for_HDL <- function(sumstats, dict) {
  
  sumstats <- sumstats[sumstats$RSID  %in% dict$rsid, ]
  sumstats %>% rename()
  
  
}