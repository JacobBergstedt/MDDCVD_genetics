
fit_MVMR <- function(df_MDD, keys, mediators, CVD, covariate_group) {

  read_mediator <- function(x, keys) {
    
    keys_mediator <- keys %>% filter(Trait == x)
    
    read_tsv(paste0("./Data/Sumstats/", x, "_sumstats")) %>% 
      mutate(rsid = RSID, pval = P) %>% 
      bind_cols(select(keys_mediator, units, prevalence, phenotype = Trait))
    
   
  }
  
  keys_CVD <- keys %>% filter(Trait == CVD)
  
  df_outcome <- read_tsv(paste0("./Data/Sumstats/", CVD,"_sumstats")) %>% 
    bind_cols(select(keys_CVD, units, prevalence, phenotype = Trait))
  
  df_mediators <- map(mediators, read_mediator, keys = keys)
  
  snps_in_common <- reduce(c(list(df_MDD, df_outcome), df_mediators), 
                           ~ inner_join(select(.x, RSID), select(., RSID)))
  
  df_MDD <- df_MDD %>% filter(RSID %in% snps_in_common$RSID)
  df_mediators <- map(df_mediators, ~ filter(., RSID %in% snps_in_common$RSID))
  
  
  MDD_CLUMP <- ieugwasr::ld_clump(dat = df_MDD,
                                  clump_kb = 5000,
                                  clump_r2 = 0.001,
                                  clump_p = 5e-8,
                                  plink_bin = plinkbinr::get_plink_exe(),
                                  bfile = "./Data/REF/g1000_eur/g1000_eur")
  
  
  df_mediators_CLUMP <- map(df_mediators, ~ ieugwasr::ld_clump(.,
                                                               clump_kb = 5000,
                                                               clump_r2 = 0.001,
                                                               clump_p = 5e-8,
                                                               plink_bin = plinkbinr::get_plink_exe(),
                                                               bfile = "./Data/REF/g1000_eur/g1000_eur"))
  
  iv_set <- bind_rows(c(list(MDD_CLUMP), df_mediators_CLUMP)) %>% 
    arrange(CHR, POS) %>% 
    select(-id) %>% 
    group_by(RSID) %>% 
    top_n(1, -log10(P)) %>% 
    ungroup()
  
  iv_set_CLUMP <- ieugwasr::ld_clump(dat = iv_set,
                                     clump_kb = 5000,
                                     clump_r2 = 0.001,
                                     clump_p = 5e-8,
                                     plink_bin = plinkbinr::get_plink_exe(),
                                     bfile = "./Data/REF/g1000_eur/g1000_eur")
  
  
  iv_set_CLUMP_formatted <- format_data_TSD(iv_set_CLUMP, type = "exposure")
  
  
  # fix outcome
  df_outcome$phenotype <- CVD
  df_outcome_formatted <- format_data_TSD(df_outcome, type = "outcome", snps = iv_set_CLUMP_formatted$SNP)
  
  df_harmonised <- harmonise_data(exposure_dat = iv_set_CLUMP_formatted,
                                  outcome_dat = df_outcome_formatted) %>% 
    filter(mr_keep) %>% 
    arrange(pval.exposure) %>% 
    steiger_filtering() %>%
    filter(!(!steiger_dir & steiger_pval < 0.10)) %>% 
    as_tibble()
  
  df <- tibble(RSID = df_harmonised$SNP)
  
  mediator_names <- str_extract(mediators, "^[:alpha:]+")
  df_list <- c(list(MDD = df_MDD), set_names(df_mediators, mediator_names))
  df_list <- map2(df_list, names(df_list), ~ select(.x, RSID, "{.y}_B" := B, "{.y}_SE" := SE))
  df <- reduce(df_list, .f = ~ inner_join(.x, .y), .init = df)
  df <- inner_join(df, select(df_outcome, RSID, "{CVD}_B" := B, "{CVD}_SE" := SE))
  
  formatted_df <- format_mvmr(BXGs = select(df, all_of(paste0(c("MDD", mediator_names), "_B"))),
                              seBXGs = select(df, all_of(paste0(c("MDD", mediator_names), "_SE"))),
                              BYG = select(df, all_of(paste0(CVD, "_B"))),
                              seBYG = select(df, all_of(paste0(CVD, "_SE"))),
                              RSID = df[, c("RSID")])
  
  res <- as.data.frame(ivw_mvmr(formatted_df)) %>% 
    rownames_to_column("ID") %>% 
    inner_join(tibble(ID = paste0("exposure", 1:(length(mediator_names) + 1)),
                      Exposure = c("MDD", mediator_names))) %>% 
    mutate(outcome = CVD, 
           nsnp = nrow(df), 
           model = paste0(CVD, " ~ MDD + ", paste0(mediator_names, collapse = " + ")), 
           Instrument_vars = paste0("MDD, ", paste0(mediator_names, collapse = ", ")),
           Covariate = covariate_group)
  
  formatted_df2 <- format_mvmr(BXGs = select(df, MDD_B),
                               seBXGs = select(df, MDD_SE),
                               BYG = select(df, paste0(CVD, "_B")),
                               seBYG = select(df, paste0(CVD, "_SE")),
                               RSID = df[, c("RSID")])
  
  res2 <- as.data.frame(ivw_mvmr(formatted_df2)) %>% 
    rownames_to_column("ID") %>% 
    mutate(outcome = CVD, nsnp = nrow(df), 
           Exposure = "MDD", 
           model = paste0(CVD, " ~ MDD"),
           Instrument_vars = paste0("MDD, ", paste0(mediator_names, collapse = ", ")),
           Covariate = "None")
  
  bind_rows(res, res2)
  
}


fit_MVMR_individual_mediator <- function(df_MDD, keys, mediator, CVD) {
  
  keys_mediator <- keys %>% filter(Trait == mediator)
  keys_CVD <- keys %>% filter(Trait == CVD)
  
  df_mediator <- read_tsv(paste0("./Data/Sumstats/", mediator ,"_sumstats"))
  df_outcome <- read_tsv(paste0("./Data/Sumstats/", CVD,"_sumstats")) %>% 
    bind_cols(select(keys_CVD, units, prevalence, phenotype = Trait))
  
  snps_in_common <- select(df_MDD, RSID) %>% 
    inner_join(select(df_mediator, RSID)) %>% 
    inner_join(select(df_outcome, RSID)) %>% 
    pull(RSID)
  
  df_MDD <- df_MDD %>% filter(RSID %in% snps_in_common)
  
  df_mediator <- df_mediator %>% 
    filter(RSID %in% snps_in_common) %>% 
    bind_cols(select(keys_mediator, units, prevalence, phenotype = Trait)) %>% 
    mutate(rsid = RSID, pval = P)
  
  MDD_CLUMP <- ieugwasr::ld_clump(dat = df_MDD,
                                  clump_kb = 5000,
                                  clump_r2 = 0.001,
                                  clump_p = 5e-8,
                                  plink_bin = plinkbinr::get_plink_exe(),
                                  bfile = "./Data/REF/g1000_eur/g1000_eur")
  
  
  
  mediator_CLUMP <- ieugwasr::ld_clump(dat = df_mediator,
                                       clump_kb = 5000,
                                       clump_r2 = 0.001,
                                       clump_p = 5e-8,
                                       plink_bin = plinkbinr::get_plink_exe(),
                                       bfile = "./Data/REF/g1000_eur/g1000_eur")
  
  iv_set <- bind_rows(MDD_CLUMP, mediator_CLUMP) %>% 
    arrange(CHR, POS) %>% 
    select(-id) %>% 
    group_by(rsid) %>% 
    top_n(1, -log10(P)) %>% 
    ungroup()
  
  iv_set_CLUMP <- ieugwasr::ld_clump(dat = iv_set,
                                     clump_kb = 5000,
                                     clump_r2 = 0.001,
                                     clump_p = 5e-8,
                                     plink_bin = plinkbinr::get_plink_exe(),
                                     bfile = "./Data/REF/g1000_eur/g1000_eur")
  
  
  iv_set_CLUMP_formatted <- format_data_TSD(iv_set_CLUMP, type = "exposure")
  
  
  # Fix outcome
  df_outcome$phenotype <- CVD
  df_outcome_formatted <- format_data_TSD(df_outcome, type = "outcome", snps = iv_set_CLUMP_formatted$SNP)
  
  df_harmonised <- harmonise_data(exposure_dat = iv_set_CLUMP_formatted, outcome_dat = df_outcome_formatted) %>% 
    filter(mr_keep) %>% 
    arrange(pval.exposure) %>% 
    steiger_filtering() %>%
    filter(!(!steiger_dir & steiger_pval < 0.10)) %>% 
    as_tibble()
  
  combine <- tibble(RSID = df_harmonised$SNP) %>% 
    inner_join(select(df_MDD, RSID, MDD_B = B, MDD_SE = SE))
  
  outcome_sub <- df_outcome[, c("RSID","B","SE")]
  colnames(outcome_sub) <- c("RSID", paste0(CVD,"_B"), paste0(CVD, "_SE"))
  combine <- inner_join(combine, outcome_sub, by = "RSID")
  mediator_sub <- df_mediator[, c("RSID","B","SE")]
  colnames(mediator_sub) <- c("RSID", paste0(mediator,"_B"), paste0(mediator,"_SE"))
  combine <- inner_join(combine, mediator_sub, by = "RSID")
  
  formatted_df <- format_mvmr(BXGs = combine %>% select(all_of(c("MDD_B", paste0(mediator, "_B")))),
                              seBXGs = combine %>% select(all_of(c("MDD_SE", paste0(mediator, "_SE")))),
                              BYG = combine %>% select(all_of(c(paste(CVD, "_B", sep = "")))),
                              seBYG = combine %>% select(all_of(c(paste(CVD, "_SE", sep = "")))),
                              RSID = combine[, "RSID"])
  
  res <- as.data.frame(ivw_mvmr(formatted_df)) %>% 
    rownames_to_column("ID") %>% 
    inner_join(tibble(ID = c("exposure1", "exposure2"), Exposure = c("MDD", mediator))) %>% 
    mutate(outcome = CVD, 
           nsnp = nrow(combine), 
           model = paste0(CVD, " ~ MDD + ", mediator), 
           Instrument_vars = paste0("MDD, ", mediator),
           Covariate = mediator)
  
  formatted_df2 <- format_mvmr(BXGs = combine %>% select("MDD_B"),
                               seBXGs = combine %>% select("MDD_SE"),
                               BYG = combine %>% select(all_of(c(paste(CVD, "_B", sep = "")))),
                               seBYG = combine %>% select(all_of(c(paste(CVD, "_SE", sep = "")))),
                               RSID = combine[, c("RSID")])
  
  res2 <- as.data.frame(ivw_mvmr(formatted_df2)) %>% 
    rownames_to_column("ID") %>% 
    mutate(outcome = CVD, 
           nsnp = nrow(combine), 
           model = paste0(CVD, " ~ MDD"), 
           Exposure = "MDD", 
           Instrument_vars = paste0("MDD, ", mediator),
           Covariate = "None")
  
  bind_rows(res, res2)
  
  
}

fit_MVMR_only_MDD <- function(df_MDD, keys, CVD) {
  
  
  keys_CVD <- keys %>% filter(Trait == CVD)
  
  
  df_outcome <- read_tsv(paste0("./Data/Sumstats/", CVD,"_sumstats")) %>% 
    bind_cols(select(keys_CVD, units, prevalence, phenotype = Trait))
  
  snps_in_common <- select(df_MDD, RSID) %>% 
    inner_join(select(df_outcome, RSID)) %>% 
    pull(RSID)
  
  df_MDD <- df_MDD %>% filter(RSID %in% snps_in_common)
  
  
  MDD_CLUMP <- ieugwasr::ld_clump(dat = df_MDD,
                                  clump_kb = 5000,
                                  clump_r2 = 0.001,
                                  clump_p = 5e-8,
                                  plink_bin = plinkbinr::get_plink_exe(),
                                  bfile = "./Data/REF/g1000_eur/g1000_eur")
  
  
  iv_set_CLUMP <- MDD_CLUMP %>% 
    arrange(CHR, POS) %>% 
    select(-id) %>% 
    group_by(rsid) %>% 
    top_n(1, -log10(P)) %>% 
    ungroup()
  
  iv_set_CLUMP_formatted <- format_data_TSD(iv_set_CLUMP, type = "exposure")
  
  
  # Fix outcome
  df_outcome$phenotype <- CVD
  df_outcome_formatted <- format_data_TSD(df_outcome, type = "outcome", snps = iv_set_CLUMP_formatted$SNP)
  
  df_harmonised <- harmonise_data(exposure_dat = iv_set_CLUMP_formatted, outcome_dat = df_outcome_formatted) %>% 
    filter(mr_keep) %>% 
    arrange(pval.exposure) %>% 
    steiger_filtering() %>%
    filter(!(!steiger_dir & steiger_pval < 0.10)) %>% 
    as_tibble()
  
  combine <- tibble(RSID = df_harmonised$SNP) %>% 
    inner_join(select(df_MDD, RSID, MDD_B = B, MDD_SE = SE))
  
  outcome_sub <- df_outcome[, c("RSID","B","SE")]
  colnames(outcome_sub) <- c("RSID", paste0(CVD,"_B"), paste0(CVD, "_SE"))
  combine <- inner_join(combine, outcome_sub, by = "RSID")
  
  
  formatted_df <- format_mvmr(BXGs = select(combine, "MDD_B"),
                              seBXGs = select(combine, "MDD_SE"),
                              BYG = select(combine, all_of(paste0(CVD, "_B"))),
                              seBYG = select(combine, all_of(paste0(CVD, "_SE"))),
                              RSID = combine[, "RSID"])
  
  res <- as.data.frame(ivw_mvmr(formatted_df)) %>% 
    rownames_to_column("ID") %>% 
    mutate(outcome = CVD, 
           Exposure = "MDD",
           nsnp = nrow(combine), 
           model = paste0(CVD, " ~ MDD"), 
           Instrument_vars = "MDD",
           Covariate = "None")
  
  res
  
  
  
  
}
