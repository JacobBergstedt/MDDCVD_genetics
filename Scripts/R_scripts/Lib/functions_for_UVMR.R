
compareNA <- function(v1, v2) {
  
  same <- (v1 == v2) | (is.na(v1) & is.na(v2))
  same[is.na(same)] <- FALSE
  same
  
}


Isq <- function(y,s){
  
  k           <-  length(y)
  w           <-  1/s^2
  sum.w  <- sum(w)
  mu.hat     <- sum(y*w)/sum.w  
  Q          <- sum(w*(y-mu.hat)^2)
  Isq        <- (Q - (k-1))/Q
  Isq        <-  max(0, Isq)
  return(Isq)
  
}

format_data_TSD <- function(df, type, snps = NULL) {
  
  
  out <- TwoSampleMR::format_data(df,
                                  type = type, 
                                  header = TRUE, 
                                  snp_col = "RSID", 
                                  snps = snps,
                                  beta_col = "B",
                                  se_col = "SE", 
                                  id_col = "phenotype",
                                  effect_allele_col = "A1",
                                  other_allele_col = "A2",
                                  pval_col = "P",
                                  min_pval = 1e-300, 
                                  chr_col = "CHR",
                                  pos_col = "POS",
                                  phenotype_col = "phenotype",
                                  ncontrol_col = "N_CON",
                                  ncase_col = "N_CAS",
                                  samplesize_col = "N", 
                                  units_col = "units",
                                  eaf_col = "EAF", 
                                  log_pval = FALSE)
  
  # prev_tib <- df %>% 
  #   select(phenotype, prevalence) %>% 
  #   distinct(phenotype, .keep_all = TRUE)
  # 
  # if (type == "exposure") {
  #   
  #   left_join(out, 
  #             rename(prev_tib, exposure = phenotype, prevalence.exposure = prevalence))
  #   
  # } else if (type == "outcome") {
  #   
  #   left_join(out, 
  #             rename(prev_tib, outcome = phenotype, prevalence.outcome = prevalence))
  #   
  # }
  
  out
  
}


prepare_exposure <- function(exposure, keys, r2_clumping = 0.001) {
  
  
  keys_exposure <- keys %>% filter(Trait == exposure)
  read_exposure <- read_tsv(paste0("Data/Sumstats/", exposure, "_sumstats")) %>% 
    mutate(units = keys_exposure$units, phenotype = exposure)
  
  exposure_format <- format_data_TSD(read_exposure, type = "exposure")
  exposure_format$prevalence.exposure <- keys_exposure$prevalence
  exposure_format$rsid <- exposure_format$SNP
  exposure_format$pval <- exposure_format$pval.exposure
  exposure_format$N <- max(read_exposure$N, na.rm = TRUE)
  exposure_format_clumped <- ieugwasr::ld_clump(dat = exposure_format,
                                                clump_kb = 5000,
                                                clump_r2 = r2_clumping,
                                                clump_p = 1e-3, 
                                                plink_bin = plinkbinr::get_plink_exe(),
                                                bfile = "Data/REF/g1000_eur/g1000_eur")
  exposure_format_clumped
  
}

UVMR <- function(exposure_format_clumped, keys, outcome, steiger_filter = TRUE, p_clumping = 5e-8, r2_clumping = 0.001){
  
  exposure <- exposure_format_clumped$exposure[1]
  keys_outcome <- keys %>% filter(Trait == outcome)
  read_outcome <- read_tsv(paste0("Data/Sumstats/", outcome, "_sumstats"))
  read_outcome$phenotype <- outcome
  read_outcome$units <- keys_outcome$units
  outcome_format <- format_data_TSD(read_outcome, type = "outcome", snps = exposure_format_clumped$SNP)
  
  outcome_format$prevalence.outcome <- keys_outcome$prevalence
  df_harmonised_nonthresh <- TwoSampleMR::harmonise_data(exposure_dat = exposure_format_clumped, outcome_dat = outcome_format)
  df_harmonised_nonthresh$id.exposure <- exposure
  df_harmonised_nonthresh$id.outcome <- outcome
  df_harmonised_nonthresh <- df_harmonised_nonthresh %>% 
    filter(mr_keep) %>% 
    arrange(pval.exposure) %>% 
    steiger_filtering()
   
  if (steiger_filter) {
    
    df_harmonised_nonthresh <- df_harmonised_nonthresh %>%
      filter(!(!steiger_dir & steiger_pval < 0.10))
    
  } 
  
  df_harmonised <- df_harmonised_nonthresh %>% filter(pval.exposure < p_clumping)
  if (nrow(df_harmonised) < 10) df_harmonised <- df_harmonised_nonthresh %>% slice(1:10)
  main_result <- TwoSampleMR::mr(df_harmonised)
  main_result$Isq_unweighted <- Isq(df_harmonised$beta.exposure, df_harmonised$se.exposure)
  main_result$Isq_weighted <- Isq(df_harmonised$beta.exposure / df_harmonised$se.outcome, df_harmonised$se.exposure / df_harmonised$se.outcome) #weighted
  N <-  max(exposure_format_clumped$N)
  
  if (compareNA(df_harmonised$units.exposure[1], "log odds")) {
    r <- TwoSampleMR::get_r_from_lor(df_harmonised$beta.exposure, 
                                     df_harmonised$eaf.exposure, 
                                     ncase = df_harmonised$ncase.exposure, 
                                     ncontrol = df_harmonised$ncontrol.exposure, 
                                     prevalence = df_harmonised$prevalence.exposure)
  } else r <- TwoSampleMR::get_r_from_bsen(b = df_harmonised$beta.exposure, se = df_harmonised$se.exposure, n = N)
  
  r2  <-  sum(r^2)
  k  <-  max(main_result$nsnp)
  main_result$f <- ((N-k-1)*r2)/(k*(1-r2))
  main_result$r2 <- r2
  main_result$Steiger_filtered <- steiger_filter
  pleiotropy <- TwoSampleMR::mr_pleiotropy_test(df_harmonised)
  heterogeneity <- TwoSampleMR::mr_heterogeneity(df_harmonised, method_list = "mr_ivw")
  leave_one_out <- TwoSampleMR::mr_leaveoneout(df_harmonised)
  
  list(pleiotropy  = pleiotropy, 
       heterogeneity = heterogeneity, leave_one_out = leave_one_out, res = main_result, Instruments = df_harmonised)
  
}




