
fix_keys <- function() {
  
  wd <- "/cluster/p/p33/cluster/users/jacobb/MDDCVD_sumstats"
  keys <- read_csv("./Data/sumstats_keys.csv") %>% 
    mutate(N = ifelse(is.na(N), N_cases + N_controls, N))
  
  keys <- mutate(keys, 
                 Input_path = file.path(Prefix, paste0(File_name, ".sumstats.gz")),
                 Output_path = file.path(wd, "Data", "Sumstats", paste0(Trait, "_sumstats")))
  
  write_csv(keys, "./Data/sumstats_keys.csv")  
  
}

add_UKB_EAF <- function(df, ukb_maf) {
  
  df %>%
    inner_join(ukb_maf, by = c("SNP" = "RSID")) %>% 
    filter((A1 == A1_UKB | A1 == A2_UKB), (A2 == A1_UKB | A2 == A2_UKB)) %>% 
    mutate(EAF_UKB = ifelse(A1 == Minor_allele_UKB, MAF_UKB, 1 - MAF_UKB))
  
  
}


process_sumstats <- function(traits, keys, n_cores) {
  
  
  
  process_individual_sumstat <- function(i, keys) {
    
    
    sumstats <- read_tsv(keys$Input_path[i], col_types = c(CHR = "c")) %>% 
      filter(!CHR %in% c("X", "Y"))
    
    cols <- colnames(sumstats)
    
    
    if (grepl("^Dep_symptoms", keys$Trait[i])) {
      
      
      sumstats$P <- exp(sumstats$P)
      
    }
      
      
    if (!"CaseN" %in% cols) {
      
      sumstats <- sumstats %>% mutate(CaseN = keys$N_cases[i])
      
    }
    
    if (!"ControlN" %in% cols) {
      
      sumstats <- sumstats %>% mutate(ControlN = keys$N_controls[i])
      
    }
    
    if ((!"N" %in% cols) & ("CaseN" %in% cols)) {
      
      sumstats <- sumstats %>% mutate(N = CaseN + ControlN)
      
    } else if (!"N" %in% cols) {
      
      sumstats <- sumstats %>% mutate(N = keys$N[i])
      
    }
    
    sumstats <- rename(sumstats, N_CON = ControlN, N_CAS = CaseN)
    sumstats <- mutate(sumstats, N_EFF = 4 / (1 / N_CON + 1 / N_CAS))
    if ("EffectAllele" %in% cols) sumstats <- rename(sumstats, A1 = EffectAllele, A2 = OtherAllele)
    if ("BETA" %in% cols) sumstats <- rename(sumstats, B = BETA)
    
    out <- keys$Output_path[i]
    vroom_write(sumstats, file = out)
    
  }
  
  # fix_keys()
  keys <- keys %>% filter(Trait %in% traits)
  n_sumstats <- nrow(keys)
  mclapply(1:n_sumstats, process_individual_sumstat, keys = keys, mc.cores = n_cores)
  
  
}
