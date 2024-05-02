setwd("/cluster/p/p33/cluster/users/jacobb/MDDCVD_sumstats/")
.libPaths("/cluster/p/p33/cluster/users/jacobb/R/x86_64-pc-linux-gnu-library/4.1")

rm(list = ls())
library(TwoSampleMR)
library(dplyr)
library(data.table)
library(parallel)
source("./Scripts/R_scripts/Lib/functions_for_UVMR.R")
source("Scripts/R_scripts/Lib/functions_for_logistics.R")

#n_cores <- 8
CVD_factor_set <- get_MDD_var_list_updated_CVD()
MDD_set <- c("MDD_PGC3_V2","MDD_PGC3_noUKB_v2")
#grid <- expand.grid(CVD_var = CVD_factor_set, MDD_var = MDD_set)
#mclapply(1:48, function(i) UVMR(exposure = grid$CVD_var[i], outcome = grid$MDD_var[i], md_version = grid$MDD_var[i]), mc.cores = n_cores)





for (i in MDD_set) {
  for (j in CVD_factor_set) {
    UVMR(exposure = j,outcome = i, md_version = i)
  }
}


exposure <- c("Overall_physical_activity_Doherty_2018_V2","PAD_2021_V2","CM_2021_V2","IL6")
outcome <- "MDD_PGC3_V2"
md_version <- "MDD_PGC3_V2"

for (i in exposure) {
  read_exposure <- data.table::fread(paste0("/cluster/p/p33/cluster/users/jacobb/MDDCVD_sumstats/Data/Sumstats/",i,"_sumstats"), header = TRUE, nThread = 10)
  read_exposure$phenotype <- i
  exposure_format <- TwoSampleMR::format_data(read_exposure, type = "exposure",
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
                                        clump_p = 1,
                                        plink_bin = plinkbinr::get_plink_exe(),
                                        bfile = "Data/REF/g1000_eur/g1000_eur")
  read_outcome <- data.table::fread(paste0("/cluster/p/p33/cluster/users/jacobb/MDDCVD_sumstats/Data/Sumstats/",outcome,"_sumstats"), header = TRUE, nThread = 10)
  read_outcome$phenotype <- outcome
  outcome_format <- TwoSampleMR::format_data(read_outcome, type = "outcome",
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
  df_harmonised$id.exposure <- i
  df_harmonised$id.outcome <- outcome
  steiger <- TwoSampleMR::steiger_filtering(df_harmonised)
  df_harmonised <- subset(steiger, steiger$steiger_dir==TRUE & mr_keep==TRUE)
  df_harmonised <- df_harmonised[order(df_harmonised$p, decreasing=FALSE),]
  df_harmonised <- df_harmonised[1:10,]
  main_result <- TwoSampleMR::mr(df_harmonised) #, method_list = c("mr_ivw","mr_egger_regression","mr_weighted_median","mr_simple_mode","mr_weighted_mode","mr_wald_ratio","mr_raps")
  
  main_result$Isq_unweighted <- Isq(df_harmonised$beta.exposure,df_harmonised$se.exposure) #unweighted
  main_result$Isq_weighted <- Isq((df_harmonised$beta.exposure/df_harmonised$se.outcome),(df_harmonised$se.exposure/df_harmonised$se.outcome)) #weighted
  N=max(read_exposure$N)
  r = TwoSampleMR::get_r_from_bsen(b = df_harmonised$beta.exposure,se = df_harmonised$se.exposure,n = N)
  r2 = sum(r^2)
  k = max(main_result$nsnp)
  main_result$f <- ((N-k-1)*r2)/(k*(1-r2))
  main_result$nsnp_after_clump <- nrow(exposure_format)
  pleiotropy <- TwoSampleMR::mr_pleiotropy_test(df_harmonised)
  heterogeneity <- TwoSampleMR::mr_heterogeneity(df_harmonised,method_list = c("mr_ivw"))
  leave_one_out <- TwoSampleMR::mr_leaveoneout(df_harmonised)
  out.fname = paste0("/cluster/p/p33/cluster/users/jacobb/MDDCVD_sumstats/Results/MR/",md_version,"/UVMR.",md_version)
  write.table(as.data.frame(main_result), paste0(out.fname,".main.txt"), row.names=F,quote=F,col.names=T, append = T, sep = "\t")
  write.table(as.data.frame(pleiotropy), paste0(out.fname,".pleiotropy.txt"), row.names=F,quote=F,col.names=T, append = T, sep = "\t")
  write.table(as.data.frame(heterogeneity), paste0(out.fname,".heterogeneity.txt"), row.names=F,quote=F,col.names=T, append = T, sep = "\t")
  write.table(as.data.frame(leave_one_out), paste0(out.fname,".mr.leaveoneout.txt"), row.names=F,quote=F,col.names=T, append = T, sep = "\t")
}