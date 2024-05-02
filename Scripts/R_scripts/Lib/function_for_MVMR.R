MVMR <- function(MDD_version,Mediator_list,CVD_version){
MDD <- fread(paste0("./Data/Sumstats/",MDD_version,"_sumstats"),header = T)
MDD$rsid <- MDD$RSID
MDD$p <- MDD$P

MDD_format <- ieugwasr::ld_clump(dat = MDD,
                                 clump_kb = 5000,
                                 clump_r2 = 0.001,
                                 clump_p = 5e-8,
                                 plink_bin = plinkbinr::get_plink_exe(),
                                 bfile = "./Data/REF/g1000_eur/g1000_eur")

MDD_sub <- MDD[,c("RSID","B","SE")]
names(MDD_sub) <- c("RSID","MDD_B","MDD_SE")

###mediator###

mediator <- Mediator_list

###cvd###

CVD <- paste0(CVD_version,"_sumstats")
for (j in 1:length(mediator)) {
  iv_set <- data.frame(MDD_format$RSID,MDD_format$P,MDD_format$B,MDD_format$SE)
  names(iv_set) <- c("rsid","p","B","SE")
  for (k in 1:length(mediator[[j]])) {
    read_mediator <- fread(paste0("./Data/Sumstats/",mediator[[j]][k]), header = TRUE)
    
    read_mediator$rsid <- read_mediator$RSID
    read_mediator$p <- read_mediator$P
    
    read_mediator <- ieugwasr::ld_clump(dat = read_mediator,
                                        clump_kb = 5000,
                                        clump_r2 = 0.001,
                                        clump_p = 5e-8,
                                        plink_bin = plinkbinr::get_plink_exe(),
                                        bfile = "./Data/REF/g1000_eur/g1000_eur")
    temp <- data.frame(read_mediator$RSID,read_mediator$P,read_mediator$B,read_mediator$SE)
    names(temp) <- c("rsid","p","B","SE")
    iv_set <- rbind(iv_set, temp)}
  
  iv_set <- iv_set[order(iv_set$p, decreasing=FALSE), ]
  iv_set <- subset(iv_set, !duplicated(rsid))
  iv_set <- ieugwasr::ld_clump(dat = iv_set,
                               clump_kb = 5000,
                               clump_r2 = 0.001,
                               clump_p = 5e-8,
                               plink_bin = plinkbinr::get_plink_exe(),
                               bfile = "./Data/REF/g1000_eur/g1000_eur")
  iv_set <- inner_join(iv_set,MDD[,c("rsid","A1","A2","N","CHR","POS","N_CAS","N_CON")],by="rsid")
  
  iv_set_format <- format_data(iv_set, type = "exposure",
                               header = TRUE,
                               snp_col = "rsid",
                               snps = NULL,
                               beta_col = "B",
                               se_col = "SE",
                               effect_allele_col = "A1",
                               other_allele_col = "A2",
                               chr_col = "CHR",
                               pos_col = "POS",
                               ncase_col = "N_CAS",
                               samplesize_col = "N",
                               pval_col = "p",
                               min_pval = 1e-200,
                               log_pval = FALSE)
  
  
  for (i in CVD) {
    read_outcome <- fread(paste0("./Data/Sumstats/",i), header = TRUE)
    
    outcome_format <- format_data(read_outcome, type = "outcome",
                                  snps = iv_set_format$SNP,
                                  header = TRUE,
                                  snp_col = "RSID",
                                  beta_col = "B",
                                  se_col = "SE",
                                  id_col = "phenotype",
                                  effect_allele_col = "A1",
                                  other_allele_col = "A2",
                                  pval_col = "PVAL", 
                                  min_pval = 1e-200,
                                  phenotype_col = "phenotype",
                                  ncase_col = "N_CAS",
                                  samplesize_col = "N",
                                  log_pval = FALSE)
    
    df_harmonised <- TwoSampleMR::harmonise_data(exposure_dat = iv_set_format,
                                                 outcome_dat = outcome_format,action = 2)
    steiger <- TwoSampleMR::steiger_filtering(df_harmonised)
    df_harmonised <- subset(steiger, steiger$steiger_dir==TRUE&mr_keep==TRUE)
    
    final_df <- filter(MDD_sub, RSID %in% df_harmonised$SNP)
    
    for (p in 1:length(mediator[[j]])) {
      read_mediator <- fread(paste0("./Data/Sumstats/",mediator[[j]][p]), header = TRUE)
      read_mediator <- read_mediator[,c("RSID","B","SE")]
      colnames(read_mediator) <- c("RSID", paste0(mediator[[j]][p],"_B"), paste0(mediator[[j]][p],"_SE"))
      final_df <- inner_join(final_df,read_mediator,by = "RSID")}
    
    
    read_outcome <- read_outcome[,c("RSID","B","SE")]
    colnames(read_outcome) <- c("RSID", paste0(i,"_B"), paste0(i,"_SE"))
    final_df_plus_outcome <- inner_join(final_df,read_outcome,by = "RSID")
    
    formatted_df <- format_mvmr(BXGs = final_df_plus_outcome %>% select(all_of(c("MDD_B",paste0(mediator[[j]],"_B")))),
                                seBXGs = final_df_plus_outcome %>% select(all_of(c("MDD_SE",paste0(mediator[[j]],"_SE")))),
                                BYG = final_df_plus_outcome %>% select(all_of(c(paste(i,"_B",sep = "")))),
                                seBYG = final_df_plus_outcome %>% select(all_of(c(paste(i,"_SE",sep = "")))),
                                RSID = final_df_plus_outcome[,c("RSID")])
    ivw <- as.data.frame(ivw_mvmr(formatted_df))
    ivw$exposure <- "MDD"
    ivw$mediator <- as.character(mediator[j])
    ivw$outcome <- i
    ivw$nsnp <- nrow(final_df_plus_outcome)
    write.table(ivw, paste0("./Results/MR/",MDD_version,"/MVMR.txt"), row.names=F,quote=F,col.names=T, append = T, sep = "\t")}
}


#################################
#################################
############none mediator

MDD <- fread(paste0("./Data/Sumstats/",MDD_version,"_sumstats"),header = T)


exposure_format <- TwoSampleMR::format_data(MDD, type = "exposure",
                                            header = TRUE,
                                            snp_col = "RSID",
                                            snps = NULL,
                                            beta_col = "B",
                                            se_col = "SE",
                                            id_col = "phenotype",
                                            effect_allele_col = "A1",
                                            other_allele_col = "A2",
                                            pval_col = "P",
                                            min_pval = 1e-200,
                                            chr_col = "CHR",
                                            pos_col = "POS",
                                            phenotype_col = "phenotype",
                                            ncase_col = "N_CAS",
                                            samplesize_col = "N",
                                            eaf_col = "EAF_1KG",
                                            log_pval = FALSE)

exposure_format$rsid <- exposure_format$SNP
exposure_format$p <- exposure_format$pval.exposure
exposure_format <- ieugwasr::ld_clump(dat = exposure_format,
                                      clump_kb = 5000,
                                      clump_r2 = 0.001,
                                      clump_p = 5e-8,
                                      plink_bin = plinkbinr::get_plink_exe(),
                                      bfile = "Data/REF/g1000_eur/g1000_eur")


CVD <- paste0(CVD_version,"_sumstats")


for (i in CVD) {
  read_outcome <- fread(paste0("./Data/Sumstats/",i), header = TRUE)
  outcome_format <- format_data(read_outcome, type = "outcome",
                                snps = exposure_format$SNP,
                                header = TRUE,
                                snp_col = "RSID",
                                beta_col = "B",
                                se_col = "SE",
                                id_col = "phenotype",
                                effect_allele_col = "A1",
                                other_allele_col = "A2",
                                pval_col = "P", 
                                min_pval = 1e-200,
                                phenotype_col = "phenotype",
                                ncase_col = "N_CAS",
                                samplesize_col = "N",
                                eaf_col = "EAF_1KG",
                                log_pval = FALSE)
  
  df_harmonised <- TwoSampleMR::harmonise_data(exposure_dat = exposure_format,
                                               outcome_dat = outcome_format,action = 2)
  steiger <- TwoSampleMR::steiger_filtering(df_harmonised)
  df_harmonised <- subset(steiger, steiger_dir==TRUE & mr_keep==TRUE)
  
  final_df <- filter(MDD_sub, RSID %in% df_harmonised$SNP)
  
  read_outcome <- read_outcome[,c("RSID","B","SE")]
  colnames(read_outcome) <- c("RSID", paste0(i,"_B"), paste0(i,"_SE"))
  final_df_plus_outcome <- inner_join(final_df,read_outcome,by = "RSID")
  
  formatted_df <- format_mvmr(BXGs = final_df_plus_outcome %>% select(all_of(c("MDD_B"))),
                              seBXGs = final_df_plus_outcome %>% select(all_of(c("MDD_SE"))),
                              BYG = final_df_plus_outcome %>% select(all_of(c(paste0(i,"_B")))),
                              seBYG = final_df_plus_outcome %>% select(all_of(c(paste0(i,"_SE")))),
                              RSID = final_df_plus_outcome[,c("RSID")])
  
  ivw <- as.data.frame(ivw_mvmr(formatted_df))
  ivw$exposure <- "MDD"
  ivw$mediator <- "none"
  ivw$outcome <- i
  ivw$nsnp <- nrow(final_df_plus_outcome)
  write.table(ivw, paste0("./Results/MR/",MDD_version,"/MVMR.txt"), row.names=F,quote=F,col.names=T, append = T, sep = "\t")}
}




