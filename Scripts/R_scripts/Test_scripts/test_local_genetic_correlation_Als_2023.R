library(LAVA)
library(tidyverse)
source("./Scripts/R_scripts/Lib/helpers.R")
source("./Scripts/R_scripts/Lib/functions_for_logistics.R")
source("./Scripts/R_scripts/Lib/functions_for_LAVA.R")

integral.p = function(integral.func, K, omega, sigma, min.iter=50000, adap.thresh=c(1e-4, 1e-6)) {

  tot.iter = min.iter * 10^(0:length(adap.thresh))
  adap.thresh = c(adap.thresh,0)  # adding dummy 0 at the end to simplify loop code
  
  p = 1; curr.iter = 0
  for (i in 1:length(tot.iter)) {
    add.iter = tot.iter[i] - curr.iter
    add.p = integral.func(K, omega, sigma, n.iter=add.iter, add.reverse = TRUE)
    p = (curr.iter*p + add.iter*add.p) / tot.iter[i]
    curr.iter = tot.iter[i]
    if (all(is.na(p)) || all(p[!is.na(p)] >= adap.thresh[i])) break
  }
  return(p)
}

bivar.cond.stats = function(draw, K, sig.xy, sig.xys, var.y) {
  m = draw[2,3] + draw[1,3] + sig.xys*(draw[1,2] + draw[1,1])
  m = m/K - sig.xy
  
  v = var.y * (draw[2,2] + 2*draw[1,2] + draw[1,1])    
  v = v / K^2
  v = ifelse(v <= 0, NA, sqrt(v))
  return(c(m,v))
}

bivariate.integral = function(K, omega, sigma, n.iter=10000, add.reverse=T) {
  
  if (!add.reverse) {
    omega.null = diag(diag(omega))
    sig.use = matrix(0,3,3); sig.use[1,1] = sigma[1,1]
    theta = matrix(0,3,3); theta[-1,-1] = omega.null*K
    
    sig.xy = sigma[1,2]
    sig.xys = sig.xy/sigma[1,1]
    var.y = sigma[2,2] - (sigma[1,2]^2)/sigma[1,1]
    
    params = apply(matrixsampling::rwishart(n.iter, K, Sigma=sig.use, Theta=theta), 3, bivar.cond.stats, K=K, sig.xy, sig.xys, var.y)   # first row is means, second is SDs
    return(conditional.norm(omega[1,2], params[1,], params[2,]))
  } else {
    p1 = bivariate.integral(K, omega, sigma, n.iter/2, add.reverse=F)
    p2 = bivariate.integral(K, omega[2:1,2:1], sigma[2:1,2:1], n.iter/2, add.reverse=F)    
    return((p1+p2)/2)
  }
}


estimate.params = function(locus, phenos=NULL) {
  if (is.null(phenos)) { phenos = locus$phenos }
  P = length(phenos); x = phenos[1:(P-1)]; y = phenos[P]
  
  # estimate
  coef = t(locus$omega[y,x] %*% solve(locus$omega[x,x])); names(coef) = x
  tau = locus$omega[y,y] - locus$omega[y,x] %*% solve(locus$omega[x,x]) %*% locus$omega[x,y]
  
  # standardise
  coef.std = NULL; for (p in x) { coef.std = c(coef.std, coef[p]*sqrt(locus$omega[p,p]/locus$omega[y,y])) }
  tau.std = tau / locus$omega[y,y]
  r2 = 1 - tau.std
  
  return(data.frame(coef=coef.std, tau=tau.std, r2=r2))
}

run.bivar = function(locus, phenos=NULL, target=NULL, adap.thresh=c(1e-4, 1e-6), p.values=T, CIs=T, param.lim=1.25) {
  if (is.null(phenos)) { 
    phenos = locus$phenos
  } else {
    phenos = as.character(phenos)
    if (any(! phenos %in% locus$phenos)) { print(paste0("Error: Invalid phenotype ID(s) provided: '",paste0(phenos[! phenos %in% locus$phenos], collapse="', '"),"'")); return(NULL) }
  }
  if (!is.null(target)) {
    target = as.character(target)
    if (length(target) > 1) { print(paste0("Error: More than one target phenotype specified")); return(NULL) }
    if (! target %in% locus$phenos) { print(paste0("Error: Invalid target phenotype specified: '", target,"'")); return(NULL) }
    if (! target %in% phenos) { phenos = c(phenos,target) }		# append target to phenos if not already present
  }
  P = length(phenos)
  if (P < 2) { print(paste0("Error: Less than 2 phenotypes provided for bivariate analysis in locus: ",locus$id)); return(NULL) }
  
  if (is.null(target)) {
    pairs = t(combn(phenos,2))	# all unique phenotype pairs
  } else {
    pairs = cbind(phenos[!phenos%in%target], target)
  }
  bivar = list(); params = c("coef","r2"); ci.params = c("rho.lower","rho.upper","r2.lower","r2.upper")
  bivar = data.frame(matrix(NA, nrow(pairs), length(params)+7)); colnames(bivar) = c("phen1","phen2",params,ci.params,"p"); bivar[,c("phen1","phen2")] = pairs
  
  for (i in 1:nrow(pairs)) {
    est.params = estimate.params(locus, phenos=pairs[i,])	# estimate params
    for (p in params) { bivar[[p]][i] = signif(est.params[[p]], 6) } # store
    
    # confidence intervals
    if (CIs) {
      ci = ci.bivariate(K = locus$K, omega = locus$omega[pairs[i,],pairs[i,]], sigma = locus$sigma[pairs[i,],pairs[i,]])
      for (p in ci.params) { bivar[[p]][i] = ci[[p]] }
    }
    # p-values
    if (p.values) { bivar$p[i] = signif(integral.p(bivariate.integral, K = locus$K, omega = locus$omega[pairs[i,],pairs[i,]], sigma = locus$sigma[pairs[i,],pairs[i,]], adap.thresh=adap.thresh), 6) }
  }
  # filter any estimates that are too far out of bounds
  bivar = filter.params(data = bivar, locus.id = locus$id, params = c(params, ci.params, "p"), param.lim = param.lim)	# first param in data set must be gamma/rho; params argument just needs to list those that will be set to NA if rho is too far out of bounds
  # cap out of bounds values (CI's are already capped)
  for (p in params) { bivar[[p]] = cap(bivar[[p]], lim=c(ifelse(p=="r2", 0, -1), 1)) }	# capping rhos at -1/1, and r2s at 0/1
  colnames(bivar)[which(colnames(bivar)=="coef")] = "rho"
  
  return(bivar[,c("phen1","phen2","rho","rho.lower","rho.upper","r2","r2.lower","r2.upper","p")])
}



analyze_locus <- function(locus_pos, input) {
  
  locus <-  process.locus(locus = locus_pos, input = input)
  if (!is_null(locus)) {
    loc_info <-  data.frame(locus = locus$id, chr = locus$chr, start = locus$start, stop = locus$stop, n.snps = locus$n.snps, n.pcs = locus$K)  
    res_uni_locus <- run.univ(locus)
    res_uni <- bind_cols(loc_info, res_uni_locus)
    res_biv <- bind_cols(loc_info, run.bivar(locus, param.lim = 10))
    list(uni = res_uni, biv = res_biv)
  } else NULL
  
}

locus <- read.loci("./Data/REF/LAVA/HLA_loci.locfile") %>% slice(nrow(.))

res <- vector(mode = "list", length = 5)
res_CVD_traits <- get_MDD_var_list_only_updated_CVD()
names(res) <- res_CVD_traits

for (trait in res_CVD_traits) {
  
  keys <- get_keys() %>%
    select(phenotype = Trait,
           cases = N_cases,
           controls =  N_controls,
           filename = Output_path) %>%
    filter(phenotype %in% c(trait, "MDD_Als_2023"))
  
  keys <- keys %>%
    mutate(filename = if_else(phenotype == "MDD_Als_2023", "Data/TMP_data/MDD_Als_2023_sumstats_HM3_incl_HLA.tsv", filename))
  
  path_input_info_file <- paste0("./Data/TMP_data/input_info_file_", trait1, ".txt")
  
  write.table(keys, file = path_input_info_file, quote = FALSE)
  
  input <- process.input(input.info.file = path_input_info_file,
                         ref.prefix = "/cluster/projects/p33/github/comorment/magma/reference/magma/g1000_eur/g1000_eur",
                         phenos = keys$phenotype,
                         sample.overlap.file = "./Data/LDSC_sample_overlap_matrix_Als_2023.txt")
  
  
  res[[trait]] <- analyze_locus(locus, input)
  
}






# "Data/TMP_data/MDD_Als_2023_sumstats_HM3_incl_HLA.tsv"