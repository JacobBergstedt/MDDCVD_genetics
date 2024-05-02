


singleTrait_likelihood <- function(theta, betX, pi1, sig1, w8s, m0, nX, bn=2^7, bins = 10, M = 1e7){
  
  M = M #1e7
  piX = abs(theta[1]);
  h2X = abs(theta[2]);
  iX = abs(theta[3]);
  
  sigX = sqrt(h2X/(piX*M))
  # Number of genotyped SNPs
  m = length(betX)
  
  Rp = iX/nX
  bX = array(betX, c(1,m)) # reshape(betXY(:,1),[1,1,m]);
  
  # Define grid for FFT
  minX = mean(bX)-(5*sd(bX));
  maxX = mean(bX)+(5*sd(bX));
  dX = (maxX-minX)/(bn-1);
  minX = minX-dX/2;
  maxX = maxX+dX/2;
  
  bXi = ceiling((bX-minX)/dX);
  bXi[bXi<1] = 1;
  bXi[bXi>bn] = bn;
  
  
  if(piX > 0.2 || piX < 1e-6 || h2X > 1 || h2X < 1e-6 || iX > 1.5 || iX < 0.5 ){
    logL = 1e10
  }else{
    min_pi1 = min(pi1)-1e-10;
    max_pi1 = max(pi1)+1e-10;
    dp = (max_pi1-min_pi1)/bins;
    pc = min_pi1 + (dp * matrix(seq(0.5,(bins-0.5),1),ncol=1))
    pix = ceiling((pi1-min_pi1)/dp);
    min_sig1 = min(sig1)-1e-10;
    max_sig1 = max(sig1)+1e-10;
    ds = (max_sig1-min_sig1)/bins;
    sc = min_sig1 + (ds * matrix(seq(0.5,(bins-0.5),1),ncol=1))
    six = ceiling((sig1-min_sig1)/ds);
    cix = pix + bins*(six-1);
    uni_cix = sort(unique(cix))
    ucix = match(uni_cix, cix)
    ixMap = match(cix, uni_cix)
    Sig1 = sc[six[ucix]];
    Pi1 = pc[pix[ucix]];
    mm = length(Sig1);
    
    Ax = aperm( array(rep((1/sigX) / Sig1, bn), c(mm,bn)), c(2,1))
    Qx = aperm( array(rep(Pi1 * piX, bn), dim = c(mm, bn)),c(2,1))
    
    
    j = 0:(bn-1)
    vi = 2*pi*(j-bn/2)/(maxX-minX-dX);
    
    Rx = array(rep(vi, mm), c(bn,mm))
    
    Lx = -m0 * ( 1 - 1 / sqrt(1 + Rx^2/Ax^2))*Qx
    Le = -(1/2)*(Rp*Rx^2)
    
    # /!\ complex numbers here!
    mf_init = -2 * log(as.complex(-1)) * ( (minX+dX/2) / (maxX-minX-dX) )*j
    mf = array(rep(mf_init, mm), dim=c(bn, mm))
    
    phi = exp(Lx+Le+mf);
    
    # In R, not possible to chose dimension for FFT, so we need a loop to do it for all rho bins
    FFT=array(NA, dim=c(bn, mm))
    for(l in 1:mm){
      FFT[,l] = fft(phi[,l])
    }
    
    FFTmod_init = (1/(maxX-minX-dX))*(as.complex(-1))^(bn*((minX+dX/2)/(maxX-minX-dX)) + j)
    FFTmod = array(rep(FFTmod_init, mm), dim=c(bn,mm))
    
    FFT0 = Re(FFT*FFTmod)
    pfE = FFT0[cbind(t(bXi), ixMap)];
    length(which(pfE<0))
    # If some, remove them & update weights to keep the correct set of SNPs
    my_w8s = w8s[pfE>0]
    pfE=pfE[pfE>0]
    
    # We use m * mean(...) to account for SNPs that may have been excluded before
    logL = -m * mean(log(pfE*my_w8s))
  }
  return(logL)
}




clump_data <- function (dat, clump_kb = 10000, clump_r2 = 0.001, clump_p1 = 1, 
          clump_p2 = 1, pop = "EUR") 
{
  

  pval_column <- "pval.exposure"
  if (!is.data.frame(dat)) {
    stop("Expecting data frame returned from format_data")
  }
  if ("pval.exposure" %in% names(dat) & "pval.outcome" %in% 
      names(dat)) {
    message("pval.exposure and pval.outcome columns present. Using pval.exposure for clumping.")
  }
  else if (!"pval.exposure" %in% names(dat) & "pval.outcome" %in% 
           names(dat)) {
    message("pval.exposure column not present, using pval.outcome column for clumping.")
    pval_column <- "pval.outcome"
  }
  else if (!"pval.exposure" %in% names(dat)) {
    message("pval.exposure not present, setting clumping p-value to 0.99 for all variants")
    dat$pval.exposure <- 0.99
  }
  else {
    pval_column <- "pval.exposure"
  }
  if (!"id.exposure" %in% names(dat)) {
    dat$id.exposure <- random_string(1)
  }
  d <- data.frame(rsid = dat$SNP, pval = dat[[pval_column]], 
                  id = dat$id.exposure)
  out <- ieugwasr::ld_clump(d, clump_kb = clump_kb, clump_r2 = clump_r2, 
                            clump_p = clump_p1, pop = pop, bfile = "/cluster/p/p33/cluster/users/jacobb/MDDCVD_sumstats/Data/REF/g1000_eur/g1000_eur", plink_bin = plinkbinr::get_plink_exe())
  
  
  keep <- paste(dat$SNP, dat$id.exposure) %in% paste(out$rsid, 
                                                     out$id)
  return(dat[keep, ])
}



calculate_SP <- function (input.df, trait.names, run_ldsc = TRUE, run_MR = TRUE, 
          saveRFiles = TRUE, hm3 = NA, ld = NA, nStep = 2, SP_single = 3, 
          SP_pair = 50, SNP_filter = 10, SNP_filter_ldsc = NA, nCores = 1, 
          M = 1e+07) 
{
  if (nStep > 2 || nStep < 1) {
    cat(print("Please choose 1 or 2 for the number of analysis steps\n"))
    stop()
  }
  if (is.na(nCores)) {
    nCores = max(1, floor((parallel::detectCores())/3))
  }
  else {
    if (nCores > parallel::detectCores()) {
      cat(print("Core number chosen is greater than cores available\n"))
      stop()
    }
  }
  if (is.na(SNP_filter_ldsc)) {
    SNP_filter_ldsc = SNP_filter
  }
  if (SNP_filter_ldsc < 1 | SNP_filter < 1) {
    cat(print("Please choose a value equal to or greater than 1 for the thinning of every nth SNP\n"))
    stop()
  }
  

  EXP = trait.names[1]
  OUT = trait.names[2]
  SP = gettingSP_ldscMR(input.df, trait.names, run_ldsc, run_MR, 
                        saveRFiles, SNP_filter_ldsc, hm3, ld)
  i_XY = as.numeric(SP[[1]])
  axy_MR = as.numeric(SP[[2]])
  ayx_MR = as.numeric(SP[[3]])
  input.df_filtered <- input.df %>% slice(seq(1, nrow(input.df), 
                                              by = SNP_filter))
  
  nX = mean(input.df_filtered$N.x)
  nY = mean(input.df_filtered$N.y)
  m0 = mean(input.df_filtered$M_LOCAL)
  bX = input.df_filtered$TSTAT.x/sqrt(nX)
  bY = input.df_filtered$TSTAT.y/sqrt(nY)
  ld = input.df_filtered$LDSC
  w8s = input.df_filtered$WEIGHT
  pi1 = input.df_filtered$PIK
  sig1 = input.df_filtered$SIGK
  sp_piX = runif(SP_single, 0, 0.01)
  sp_h2X = runif(SP_single, 0, 0.5)
  sp_iX = runif(SP_single, 0.5, 1.5)
  para = cbind(sp_piX, sp_h2X, sp_iX)
  sp_mat = matrix(unlist(para), ncol = 3, byrow = FALSE)
  colnames(sp_mat) = colnames(para)
  par.df = data.frame(par = I(apply(sp_mat, 1, as.list)))
  test.exp <- parallel::mclapply(par.df[[1]], function(x) {
    theta = unlist(x)
    test1 = optim(theta, singleTrait_likelihood, betX = bX, 
                  pi1 = pi1, sig1 = sig1, w8s = w8s, M = M, m0 = m0, 
                  nX = nX, bn = 2^7, bins = 10, method = "Nelder-Mead", 
                  control = list(maxit = 5000))
    list(mLL = test1$value, par = test1$par, conv = test1$convergence)
  }, mc.cores = nCores)
  test.exp = as.data.frame(t(matrix(unlist(test.exp), nrow = length(unlist(test.exp[1])))))
  test.out <- parallel::mclapply(par.df[[1]], function(x) {
    theta = unlist(x)
    test2 = optim(theta, singleTrait_likelihood, betX = bY, 
                  pi1 = pi1, sig1 = sig1, w8s = w8s, M = M, m0 = m0, 
                  nX = nY, bn = 2^7, bins = 10, method = "Nelder-Mead", 
                  control = list(maxit = 5000))
    list(mLL = test2$value, par = test2$par, conv = test2$convergence)
  }, mc.cores = nCores)
  test.out = as.data.frame(t(matrix(unlist(test.out), nrow = length(unlist(test.out[1])))))
  colnames(test.exp) = c("mLL", "piX", "h2X", "iX", "conv")
  colnames(test.out) = c("mLL", "piX", "h2X", "iX", "conv")
  res_exp_min = test.exp[which(test.exp$mLL == min(test.exp$mLL)), 
  ]
  res_exp = abs(res_exp_min[2:4])
  res_out_min = test.out[which(test.out$mLL == min(test.out$mLL)), 
  ]
  res_out = abs(res_out_min[2:4])
  if (saveRFiles) {
    res_singleTrait = rbind(res_exp_min, res_out_min)
    res_singleTrait = cbind(Trait = c(EXP, OUT), res_singleTrait)
    write.csv(res_singleTrait, paste0("SingleTraitAnalysis_", 
                                      EXP, "-", OUT, ".csv"), row.names = F)
  }
  pi_X = as.numeric(res_exp[1])
  pi_Y = as.numeric(res_out[1])
  h2_x = as.numeric(res_exp[2])
  h2_y = as.numeric(res_out[2])
  i_X = as.numeric(res_exp[3])
  i_Y = as.numeric(res_out[3])
  if (nStep == 2) {
    sp_tX = runif(SP_pair, 0, 0.5)
    sp_tY = runif(SP_pair, -0.5, 0.5)
    sp_h2X = max(0, h2_x - (sp_tX^2))
    sp_h2Y = max(0, h2_y - (sp_tY^2))
    sp_axy = replicate(SP_pair, (axy_MR + runif(1, -0.1, 
                                                0.1)))
    sp_ayx = replicate(SP_pair, (ayx_MR + runif(1, -0.1, 
                                                0.1)))
    sp_iXY = replicate(SP_pair, (i_XY + runif(1, -0.05, 0.05)))
    para = cbind(sp_h2X, sp_h2Y, sp_tX, sp_tY, sp_axy, sp_ayx, 
                 sp_iXY)
    sp_mat1 = matrix(unlist(para), ncol = 7, byrow = FALSE)
    colnames(sp_mat1) = c("sp_h2X", "sp_h2Y", "sp_tX", "sp_tY", 
                          "sp_axy", "sp_ayx", "sp_iXY")
    if (saveRFiles) {
      write.csv(sp_mat1, "StartingPoints.csv", row.names = F)
    }
    return(list(iX = i_X, iY = i_Y, piX = pi_X, piY = pi_Y, 
                i_XY = i_XY, axy_MR = axy_MR, ayx_MR = ayx_MR, input.df_filtered = input.df_filtered, 
                sp_mat = sp_mat1))
  }
  if (nStep == 1) {
    sp_piX = runif(SP_pair, 0, 1e-04)
    sp_piY = runif(SP_pair, 0, 1e-04)
    sp_tX = runif(SP_pair, 0, 0.5)
    sp_tY = runif(SP_pair, -0.5, 0.5)
    sp_h2X = max(0, h2_x - (sp_tX^2))
    sp_h2Y = max(0, h2_y - (sp_tY^2))
    sp_axy = replicate(SP_pair, (axy_MR + runif(1, -0.1, 
                                                0.1)))
    sp_ayx = replicate(SP_pair, (ayx_MR + runif(1, -0.1, 
                                                0.1)))
    sp_iXY = replicate(SP_pair, (i_XY + runif(1, -0.05, 0.05)))
    para = cbind(sp_piX, sp_piY, sp_h2X, sp_h2Y, sp_tX, sp_tY, 
                 sp_axy, sp_ayx, sp_iXY)
    sp_mat1 = matrix(unlist(para), ncol = 9, byrow = FALSE)
    colnames(sp_mat1) = c("sp_piX", "sp_piY", "sp_h2X", "sp_h2Y", 
                          "sp_tX", "sp_tY", "sp_axy", "sp_ayx", "sp_iXY")
    if (saveRFiles) {
      write.csv(sp_mat1, "StartingPoints.csv", row.names = F)
    }
    return(list(iX = i_X, iY = i_Y, i_XY = i_XY, axy_MR = axy_MR, 
                ayx_MR = ayx_MR, input.df_filtered = input.df_filtered, 
                sp_mat = sp_mat1))
  }
}


gettingSP_ldscMR <- function(input.df,trait.names,run_ldsc=TRUE,run_MR=TRUE,saveRFiles=TRUE,SNP_filter_ldsc,hm3,ld){
  
  
  
  
  EXP = trait.names[1]
  OUT = trait.names[2]
  nX = mean(input.df$N.x)  #get sample size for trait X
  nY = mean(input.df$N.y)  #get sample size for trait Y
  
  bX = input.df$TSTAT.x/sqrt(nX)   #get standardised beta for trait X
  bY = input.df$TSTAT.y/sqrt(nY)   #get standardised beta for trait Y
  
  X = select(input.df, RSID, CHR, A1, A2, BETA.x, SE.x, PVAL.x, TSTAT.x, N.x) %>% rename(unstdb = BETA.x, sderr = SE.x, pval = PVAL.x, TSTAT=TSTAT.x, N = N.x)
  Y = select(input.df, RSID, CHR, A1, A2, BETA.y, SE.y, PVAL.y, TSTAT.y, N.y) %>% rename(unstdb = BETA.y, sderr = SE.y, pval = PVAL.y, TSTAT=TSTAT.y, N = N.y)
  X$BETA = bX
  Y$BETA = bY
  X$SE = 1/sqrt(nX)
  Y$SE = 1/sqrt(nY)
  
  if(run_ldsc){
    # Slice, every 10th SNP for faster computation
    X %>%
      slice(seq(1, nrow(X), by=SNP_filter_ldsc)) -> X_filtered
    #nrow(X_filtered)
    Y %>%
      slice(seq(1, nrow(Y), by=SNP_filter_ldsc)) -> Y_filtered
    #nrow(Y_filtered)
    
    random_code <- paste0(LETTERS[sample.int(10, 10)], collapse = "")
    exp_file <- paste0(EXP, random_code, "_GWAS.txt")
    out_file <- paste0(OUT, random_code, "_GWAS.txt")
    write.table(X_filtered, file = exp_file, sep = "\t", quote = FALSE, row.names = FALSE, col.names=TRUE)
    write.table(Y_filtered, file = out_file, sep = "\t", quote = FALSE, row.names = FALSE, col.names=TRUE)
    
    cat("Munging exposure and outcome data FOR LDSC: \n")
    invisible(utils::capture.output(GenomicSEM::munge( exp_file,
                                                       hm3,
                                                       trait.names=paste0(EXP, random_code))))
    invisible(utils::capture.output(GenomicSEM::munge( out_file,
                                                       hm3,
                                                       trait.names=paste0(OUT, random_code))))
    
    cat("Please check the log files", paste0(EXP, "_munge.log and ", OUT, "_munge.log" ),  "to ensure that all columns were interpreted correctly and no warnings were issued for any of the summary statistics files\n")
    
    traits = c(paste0(EXP, random_code, ".sumstats.gz"), paste0(OUT, random_code, ".sumstats.gz"))
    sample.prev <- c(NA,NA) #continuous traits
    population.prev <- c(NA,NA) #continuous traits
    
    trait.names<-c(EXP, OUT)
    invisible(utils::capture.output(LDSCoutput <- GenomicSEM::ldsc(traits,
                                                                   sample.prev,
                                                                   population.prev,
                                                                   ld,
                                                                   trait.names)))
    #save(LDSCoutput, file="Pfactor.RData")
    
    cat("Please check the log file", paste0(c(traits, "ldsc.log"), collapse="_"),  "for in-depth results of the cross-trait LDSC analysis\n")
    
    # Clean up generated files, but save log file
    file.remove(paste0(paste0(EXP, random_code), c("_GWAS.txt", ".sumstats.gz")))  #need to silent the result
    file.remove(paste0(paste0(OUT, random_code), c("_GWAS.txt", ".sumstats.gz"))) #need to silent the result
    #i_X = as.numeric(LDSCoutput$I[1,1])
    #i_Y = as.numeric(LDSCoutput$I[2,2])
    i_XY = as.numeric(LDSCoutput$I[1,2])
    #h2_x_ldsc = as.numeric(LDSCoutput$S[1,1])
    #h2_y_ldsc = as.numeric(LDSCoutput$S[2,2])
  } else {i_XY = runif(1,-0.2,0.2)}
  
  
  # Running standard MR
  if(run_MR){
    MR_output = paste0(EXP,"-",OUT,"_MRresults.csv")
    # Get significant SNPs above certain Z-statistic corresponding to set p-value
    prune_X = function(zX,p_limit=1e-5){
      zX=zX
      z_limit=abs(qnorm(0.5*p_limit))
      ind_keep=which(abs(zX)>z_limit)
      ind_keep=unique(ind_keep)
      ind_keep=list(ind_keep)
      return(ind_keep)
    }
    
    # Taken from Jonathan Sulc to create bins that fit a maximum of 50k SNPs (max threshold for clumping)
    snp_bin  =  function( snp_ranks,
                          chunk_size = 50000 ){
      if (nrow( snp_ranks ) == 0) {
        return()
      }
      
      max_chr  =  snp_ranks$chr %>%
        table %>%
        cumsum %>%
        (function(x) x < chunk_size) %>%
        (function(x) names(x)[ max(which(x)) ] ) %>%
        as.numeric
      if (is.na( max_chr )) {
        max_chr = min( snp_ranks$chr )
      }
      
      bin = snp_ranks %>%
        dplyr::filter( chr <= max_chr ) %>%
        list
      return( c( bin,
                 snp_bin( snp_ranks[ snp_ranks$chr > max_chr, ],
                          chunk_size ) ) )
    }
    
    # Set threshold values
    pval=2*pnorm(-abs(5.45))
    pval1=2*pnorm(-abs(4))
    reverse_t_threshold  =  qnorm( 5e-2 )
    
    # Forward MR estimation
    mr_dataX = cbind.data.frame(SNP = X$RSID, beta = X$BETA, se = X$SE, effect_allele = X$A1, other_allele = X$A2, chr=X$CHR, Phenotype=EXP, tstat=X$TSTAT )
    mr_dataY = cbind.data.frame(SNP = Y$RSID, beta = Y$BETA, se = Y$SE, effect_allele = Y$A1, other_allele = Y$A2, chr=Y$CHR, Phenotype=OUT, tstat=Y$TSTAT )
    
    mr_ind=unlist(prune_X(mr_dataX$tstat,pval))
    print(length(mr_ind))
    if(length(mr_ind)==0){
      mr_ind=unlist(prune_X(mr_dataX$tstat,pval1))
      print(length(mr_ind))
    }
    mr_dataX = mr_dataX[mr_ind,]
    mr_dataY = mr_dataY[mr_ind,]
    
    # Remove SNPs that are more strongly associated with the outcome than the exposure
    ind_keep=which((abs(mr_dataX$beta)-abs(mr_dataY$beta))/sqrt(mr_dataX$se^2+mr_dataY$se^2) > reverse_t_threshold)
    print(length(ind_keep))
    mr_dataX = mr_dataX[ind_keep,]
    mr_dataY = mr_dataY[ind_keep,]
    
    exp_dat <- TwoSampleMR::format_data(mr_dataX, type="exposure")  #same rows as mr_dataX
    clump_bin = snp_bin(mr_dataX,50000)
    
    exp_data = c()
    for (x in 1:length(clump_bin)) {
      temp = exp_dat[exp_dat$SNP %in% clump_bin[[x]]$SNP,]
      #temp1 = TwoSampleMR::clump_data(temp)
      temp1 = clump_data(temp)
      exp_data=rbind(exp_data,temp1)
    }
    
    dups=which(duplicated(exp_data$SNP)==TRUE)
    if(length(dups)>0){
      exp_dat2 = exp_data[-dups,]
    }else{
      exp_dat2 = exp_data
    }
    
    out_dat <- TwoSampleMR::format_data(mr_dataY, type="outcome")
    out_dat2=out_dat[out_dat$SNP %in% exp_dat2$SNP,]
    
    exp_dat2=exp_dat2[order(exp_dat2$SNP),]
    out_dat2=out_dat2[order(out_dat2$SNP),]
    
    if(all(exp_dat2$`effect_allele.exposure` == out_dat2$`effect_allele.outcome`)){
      print("action=1")
      action = 1
    } else {
      print("action=2/3")
      aligned = which(exp_dat2$`effect_allele.exposure` == out_dat2$`effect_allele.outcome` &
                        exp_dat2$`other_allele.exposure` == out_dat2$`other_allele.outcome`)
      swapped = which(exp_dat2$`effect_allele.exposure` == out_dat2$`other_allele.outcome` &
                        exp_dat2$`other_allele.exposure` == out_dat2$`effect_allele.outcome`)
      exp_dat2[swapped,'beta.exposure']=exp_dat2[swapped,'beta.exposure']*-1
      exp_dat2 = exp_dat2[c(aligned,swapped),]
      out_dat2 = out_dat2[c(aligned,swapped),]
      action = 1  #made sure all strands are okay
    }
    # Harmonise the exposure and outcome data
    dat <- TwoSampleMR::harmonise_data(
      exposure_dat = exp_dat2,
      outcome_dat = out_dat2, action = action
    )
    
    # Sensitivity - Q-test
    het <- TwoSampleMR::mr_heterogeneity(dat)
    het$I2 = ((het$Q-het$Q_df)/het$Q)*100
    plei <- TwoSampleMR::mr_pleiotropy_test(dat)
    
    smaller=FALSE
    tryCatch( {res1 <- TwoSampleMR::mr(dat); print("Bigger MR list") }, error = function(e) {smaller<<-TRUE})
    if(smaller){
      print("Smaller MR list")
      res1 <- TwoSampleMR::mr(dat, method_list=c("mr_egger_regression","mr_weighted_median","mr_ivw"))
    }else{
      res1 <- TwoSampleMR::mr(dat)
    }
    
    if(nrow(res1)<2){
      axy_MR = res1[,'b']
    }else{
      axy_MR = res1[which(res1$method=="Inverse variance weighted"),'b']
    }
    
    if(saveRFiles){
      write.table(as.data.frame(res1), file=MR_output, sep = ",", append = TRUE, row.names = FALSE)
      write.table(as.data.frame(het), file=MR_output, sep = ",", append = TRUE, row.names = FALSE)
      write.table(as.data.frame(plei), file=MR_output, sep = ",", append = TRUE, row.names = FALSE)
      write.table("*", MR_output, sep = ",", quote = FALSE, col.names = FALSE, row.names = FALSE, append = TRUE)
    }
    
    # Reverse MR estimation
    rm(mr_dataX,mr_dataY,exp_dat,exp_data,exp_dat2,out_dat,out_dat2,dups,ind_keep,mr_ind,clump_bin, action, temp, temp1,dat,het,plei,smaller)
    
    # Reverse the exposure and outcome to Y - X, nothing else besides this needs to change
    mr_dataX = cbind.data.frame(SNP = Y$RSID, beta = Y$BETA, se = Y$SE, effect_allele = Y$A1, other_allele = Y$A2, chr=Y$CHR, Phenotype=OUT, tstat = Y$TSTAT )
    mr_dataY = cbind.data.frame(SNP = X$RSID, beta = X$BETA, se = X$SE, effect_allele = X$A1, other_allele = X$A2, chr=X$CHR, Phenotype=EXP, tstat = X$TSTAT )
    
    mr_ind=unlist(prune_X(mr_dataX$tstat,pval))
    print(length(mr_ind))
    if(length(mr_ind)==0){
      mr_ind=unlist(prune_X(mr_dataX$tstat,pval1))
      print(length(mr_ind))
    }
    mr_dataX = mr_dataX[mr_ind,]
    mr_dataY = mr_dataY[mr_ind,]
    
    ind_keep=which((abs(mr_dataX$beta)-abs(mr_dataY$beta))/sqrt(mr_dataX$se^2+mr_dataY$se^2) > reverse_t_threshold)
    print(length(ind_keep))
    mr_dataX = mr_dataX[ind_keep,]
    mr_dataY = mr_dataY[ind_keep,]
    
    exp_dat <- TwoSampleMR::format_data(mr_dataX, type="exposure")  #same rows as mr_dataX
    clump_bin = snp_bin(mr_dataX,50000)
    
    exp_data = c()
    for (x in 1:length(clump_bin)) {
      temp = exp_dat[exp_dat$SNP %in% clump_bin[[x]]$SNP,]
      #temp1 = TwoSampleMR::clump_data(temp)
      temp1 = clump_data(temp)
      exp_data=rbind(exp_data,temp1)
    }
    
    dups=which(duplicated(exp_data$SNP)==TRUE)
    if(length(dups)>0){
      exp_dat2 = exp_data[-dups,]
    }else{
      exp_dat2 = exp_data
    }
    
    out_dat <- TwoSampleMR::format_data(mr_dataY, type="outcome")
    out_dat2=out_dat[out_dat$SNP %in% exp_dat2$SNP,]
    
    exp_dat2=exp_dat2[order(exp_dat2$SNP),]
    out_dat2=out_dat2[order(out_dat2$SNP),]
    
    if(all(exp_dat2$`effect_allele.exposure` == out_dat2$`effect_allele.outcome`)){
      print("action=1")
      action = 1
    } else {
      print("action=2/3")
      aligned = which(exp_dat2$`effect_allele.exposure` == out_dat2$`effect_allele.outcome` &
                        exp_dat2$`other_allele.exposure` == out_dat2$`other_allele.outcome`)
      swapped = which(exp_dat2$`effect_allele.exposure` == out_dat2$`other_allele.outcome` &
                        exp_dat2$`other_allele.exposure` == out_dat2$`effect_allele.outcome`)
      exp_dat2[swapped,'beta.exposure']=exp_dat2[swapped,'beta.exposure']*-1
      exp_dat2 = exp_dat2[c(aligned,swapped),]
      out_dat2 = out_dat2[c(aligned,swapped),]
      action = 1  #made sure all strands are okay
    }
    
    # Harmonise the exposure and outcome data
    dat <- TwoSampleMR::harmonise_data(
      exposure_dat = exp_dat2,
      outcome_dat = out_dat2, action = action
    )
    
    # Sensitivity - Q-test
    het <- TwoSampleMR::mr_heterogeneity(dat)
    het$I2 = ((het$Q-het$Q_df)/het$Q)*100
    plei <- TwoSampleMR::mr_pleiotropy_test(dat)
    
    smaller=FALSE
    tryCatch( {res2 <- TwoSampleMR::mr(dat); print("Bigger MR list") }, error = function(e) {smaller<<-TRUE})
    if(smaller){
      print("Smaller MR list")
      res2 <- TwoSampleMR::mr(dat, method_list=c("mr_egger_regression","mr_weighted_median","mr_ivw"))
    }else{
      res2 <- TwoSampleMR::mr(dat)
    }
    
    if(nrow(res2)<2){
      ayx_MR = res2[,'b']
    }else{
      ayx_MR = res2[which(res2$method=="Inverse variance weighted"),'b']
    }
    
    if(saveRFiles){
      write.table(as.data.frame(res2), file=MR_output, sep = ",", append = TRUE, row.names = FALSE)
      write.table(as.data.frame(het), file=MR_output, sep = ",", append = TRUE, row.names = FALSE)
      write.table(as.data.frame(plei), file=MR_output, sep = ",", append = TRUE, row.names = FALSE)
    }
    
  } else{axy_MR = runif(1,-0.5,0.5);ayx_MR = runif(1,-0.5,0.5);}
  
  
  return(list(i_XY, axy_MR, ayx_MR))
  
}