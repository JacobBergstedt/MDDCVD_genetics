


run.pcor <- function (locus, target, phenos = NULL, adap.thresh = c(1e-04, 
                                                        1e-06), p.values = T, CIs = T, max.r2 = 0.95, param.lim = 1.25) {
  if (is.null(phenos)) {
    phenos = locus$phenos
  }
  else {
    phenos = as.character(phenos)
    if (any(!phenos %in% locus$phenos)) {
      stop(paste0("Invalid phenotype ID(s) provided: '", 
                  paste0(phenos[!phenos %in% locus$phenos], collapse = "', '"), 
                  "'"))
    }
  }
  target = as.character(target)
  if (length(target) != 2) {
    stop(paste0("Exactly two phenotype IDs must be provided as the target"))
  }
  if (!all(target %in% locus$phenos)) {
    stop(paste0("Invalid target phenotype specified: '", 
                paste0(target[!target %in% locus$phenos], collapse = "', '"), 
                "'"))
  }
  phenos = unique(c(target, phenos))
  P = length(phenos)
  if (P < 3) {
    stop(paste0("Less than 3 phenotypes provided for partial correlation analysis in locus: ", 
                locus$id))
  }
  x = which(phenos == target[1])
  y = which(phenos == target[2])
  z = which(!phenos %in% target)
  print(paste0("~ Running partial correlation for '", phenos[x], 
               "' and '", phenos[y], "', conditioned on '", paste(phenos[z], 
                                                                  collapse = "' + '"), "'"))
  eig = eigen(locus$omega[phenos, phenos][z, z])
  if (any(eig$values/(sum(eig$values)/length(z)) < 1e-04)) {
    stop("Genetic covariance matrix for phenotypes in Z is not invertible")
  }
  rm(eig)
  r2 = data.frame(x = NA, y = NA)
  for (i in 1:2) {
    if (P == 3) {
      r2[i] = run.bivar(locus, phenos[c(z, i)], p.values = F, 
                        CIs = F)$r2
    }
    else {
      r2[i] = run.multireg(locus, target = phenos[i], phenos = phenos[z], only.full.model = T, p.values = F, CIs = F, suppress.message = T)[[1]][[1]]$r2[1]
    }
  }
  out = data.frame(phen1 = phenos[x], phen2 = phenos[y], z = paste(phenos[z], 
                                                                   collapse = ";"), r2.phen1_z = r2$x, r2.phen2_z = r2$y, 
                   pcor = NA, ci.lower = NA, ci.upper = NA, p = NA)
  out$pcor = partial.cor(locus$omega[phenos, phenos], x, y, z)
  if (CIs) {
    ci = ci.pcor(K = locus$K, xy.index = c(x, y), z.index = z, 
                 omega = locus$omega[phenos, phenos], sigma = locus$sigma[phenos, phenos])
    out$ci.lower = ci[2]
    out$ci.upper = ci[3]
  }
  if (all(out[c("r2.phen1_z", "r2.phen2_z")] < max.r2) & p.values) {
    out$p = integral.p(pcov.integral, K = locus$K, omega = locus$omega[phenos, phenos], sigma = locus$sigma[phenos, phenos], adap.thresh = adap.thresh)
  }
  params = c("pcor", "ci.lower", "ci.upper")
  out = filter.params(data = out, locus.id = locus$id, params = params, param.lim = param.lim)
  out$pcor = cap(out$pcor)
  out[, !colnames(out) %in% c("phen1", "phen2", "z")] = as.data.frame(lapply(out[, !colnames(out) %in% c("phen1", "phen2", "z")], signif, 6))
  return(out)
}

pcov.cond.stats = function(draw, K, sigma, gamma, sig.xys, var.y) {
  Pw = dim(sigma)[1]-1; Pz = Pw - 1
  i.eps = 1:Pw; i.eps.z = 1:Pz; i.eps.x = Pw
  i.delta = 1:Pw + Pw; i.delta.z = 1:Pz + Pw; i.x = 2*Pw; i.y = i.x+1
  
  dtd.w = matrix(draw[i.eps,i.eps] + draw[i.eps,i.delta] + draw[i.delta,i.eps] + draw[i.delta,i.delta], ncol=Pw)
  dtd.w[Pw,Pw] = dtd.w[Pw,Pw] + 2*t(gamma$x) %*% draw[i.delta.z,i.eps.x] + gamma$fit.x
  dtd.w[-Pw,Pw] = dtd.w[-Pw,Pw] + (draw[i.delta.z,i.delta.z] + draw[i.eps.z,i.delta.z]) %*% gamma$x
  dtd.w[Pw,-Pw] = dtd.w[-Pw,Pw]
  
  omega.w = matrix(dtd.w/K - sigma[1:Pw,1:Pw], ncol=Pw)
  omega.z.inv = tryCatch(solve(omega.w[-Pw,-Pw]),error=function(x){omega.w[-Pw,-Pw]*NA})  # silently put to NA if not invertible
  b = matrix(c(-(omega.z.inv %*% omega.w[1:Pz,Pw]), 1), ncol=1)
  
  dhw.dy = gamma$dw.dz.gamma + draw[i.eps,i.delta.z] %*% gamma$y + draw[i.eps,i.y]
  dw.ew = rbind(draw[i.delta.z,i.eps], t(gamma$x) %*% draw[i.delta.z,i.eps] + draw[i.x,i.eps])
  
  M = dhw.dy/K + (dw.ew + draw[i.eps,i.eps]) %*% sig.xys - sigma[1:Pw,Pw+1]
  V = t(b) %*% dtd.w %*% b * var.y  
  
  return(c(t(b) %*% M, ifelse(V >= 0, sqrt(V), NA)))
}

# checks indices, puts xy at end
check.index = function(xy.index, z.index, P) {
  xy.index = unique(round(xy.index))
  z.index = unique(round(z.index))
  index = as.numeric(c(z.index,xy.index))
  if (any(is.na(index) | index > P | index <= 0) || length(xy.index) != 2 || length(z.index) == 0 || any(z.index %in% xy.index)) index = NULL
  return(index) 
}

estimate.pcor = function(draw, sigma) {
  i.y = dim(sigma)[1]; i.x = i.y - 1; i.z = 1:(i.x-1)
  o = draw - sigma
  
  cov = o[i.x,i.y] - t(o[i.z,i.x]) %*% solve(o[i.z,i.z]) %*% o[i.z,i.y]   
  vars = c(
    o[i.x,i.x] - t(o[i.z,i.x]) %*% solve(o[i.z,i.z]) %*% o[i.z,i.x], 
    o[i.y,i.y] - t(o[i.z,i.y]) %*% solve(o[i.z,i.z]) %*% o[i.z,i.y]
  ) 
  r = suppressWarnings(cov / sqrt(prod(vars)))
  return(r)
}

conditional.norm = function(obs, means, sds) {
  obs = abs(obs)
  prob = suppressWarnings(pnorm(obs, mean=means, sd=sds, lower.tail=F))
  prob = prob + suppressWarnings(pnorm(-obs, mean=means, sd=sds, lower.tail=T))
  return(mean(prob, na.rm=T))
}

# single indices x and y, vector of indices z
partial.cor = function(omega, x, y, z) {
  p.cov = partial.cov(omega, x, y, z)
  return(p.cov/sqrt(partial.var(omega, x, z) * partial.var(omega, y, z)))
}

# single indices x and y, vector of indices z
partial.cov = function(omega, x, y, z) {omega[x,y] - t(omega[z,x]) %*% solve(omega[z,z]) %*% omega[z,y]}

partial.var = function(omega, x, z) {omega[x,x] - t(omega[z,x]) %*% solve(omega[z,z]) %*% omega[z,x]}



ci.bivariate = function(K, omega, sigma, n.iter=10000) {  
  S = diag(sqrt(diag(omega)))
  corrs = solve(S) %*% omega %*% solve(S)
  corrs[corrs >= 1] = 0.99999; corrs[corrs <= -1] = -0.99999; diag(corrs) = 1
  omega = S %*% corrs %*% S
  
  P = dim(omega)[1]; tri = lower.tri(corrs)
  out = data.frame(pheno1=col(corrs)[tri], pheno2=row(corrs)[tri], r=corrs[tri], rho.lower=NA, rho.upper=NA, r2.lower=NA, r2.upper=NA)
  draws = tryCatch(matrixsampling::rwishart(n.iter, K, Sigma=sigma/K, Theta=omega), error=function(e){return(NULL)})
  if (!is.null(draws)) {
    if (P == 2) {
      func = function(draw, sigma) {o = draw-sigma; o[1,2]/sqrt(o[1,1]*o[2,2])}
      r = matrix(suppressWarnings(apply(draws, 3, func, sigma)), nrow=1)
    } else {
      func = function(draw, sigma) {cov2cor(draw-sigma)[lower.tri(sigma)]}
      r = suppressWarnings(apply(draws, 3, func, sigma))
    }
    # quantiles rho
    qq = round(apply(r, 1, quantile, c(0.025, 0.975), na.rm=T),5)
    qq[qq < -1] = -1; qq[qq > 1] = 1
    out$rho.lower = qq[1,]; out$rho.upper = qq[2,]
    
    # quantiles r2
    qq = round(apply(r^2, 1, quantile, c(0.025, 0.975), na.rm=T),5)
    qq[qq < 0] = 0; qq[qq > 1] = 1
    out$r2.lower = qq[1,]; out$r2.upper = qq[2,]
    if (sign(out$rho.lower)!=sign(out$rho.upper)) out$r2.lower = 0  # set r2 lower to 0 if rho CI spans 0
  }
  
  return(out)
}


### MULTILPE REG ###
# expects omega.x to be invertible
ci.multivariate = function(K, omega, sigma, n.iter=10000) {
  P = dim(omega)[1]
  S = diag(sqrt(diag(omega)))
  corrs = solve(S) %*% omega %*% solve(S)
  corrs[corrs >= 1] = 0.99999; corrs[corrs <= -1] = -0.99999; diag(corrs) = 1
  omega = S %*% corrs %*% S
  fit = omega[-P,P] %*% solve(omega[-P,-P]) %*% omega[-P,P]
  if (fit >= omega[P,P]) omega[P,P] = fit/0.99999
  # increasing omega_Y to fit if r2 > 1; setting r2 slightly below 1 in that case, since otherwise the matrixsampling::rwishart function will fail 
  
  gamma.ss = solve(corrs[-P,-P]) %*% corrs[-P,P]
  r2 = max(0, fit/omega[P,P])
  
  draws = tryCatch(matrixsampling::rwishart(n.iter, K, Sigma=sigma/K, Theta=omega), error=function(e){return(NULL)}) 
  if (!is.null(draws)) {
    est = apply(draws, 3, estimate.std, sigma)  
    qq = round(apply(est, 1, quantile, c(0.025, 0.975), na.rm=T), 5)
  } else {
    qq = matrix(NA, nrow=2, ncol=P)
  }
  qq.r2 = qq[,P]; qq.r2[qq.r2 < 0] = 0; qq.r2[qq.r2 > 1] = 1;
  
  ci = suppressWarnings(data.frame(gamma=gamma.ss, gamma.lower=qq[1,-P], gamma.upper=qq[2,-P], r2=r2, r2.lower=qq.r2[1], r2.upper=qq.r2[2]))
  return(ci)
}


estimate.std = function(draw, sigma) {
  P = dim(sigma)[1]
  o = cov2cor(draw-sigma)
  g = solve(o[-P,-P]) %*% o[-P,P]
  r2 = o[-P,P] %*% g
  return(c(g,r2))
}



### PARTIAL COR ###
# xy.index must be two unique index values, z.index any number (>0) of values not in xy.index
# expects omega.z to be invertible
# returns vector with three values: estimate (for reference), ci.low, ci.high
ci.pcor = function(K, xy.index, z.index, omega, sigma, n.iter=10000) {
  index = check.index(xy.index, z.index, dim(omega)[1]); out = rep(NA,3)
  if (is.null(index)) return(out) # just failing quietly for now
  
  omega = omega[index,index]; sigma = sigma[index,index]  
  P = dim(omega)[1]; Pw = P - 1; Pz = Pw - 1;
  S = diag(sqrt(diag(omega)))
  corrs = solve(S) %*% omega %*% solve(S)
  corrs[corrs >= 1] = 0.99999; corrs[corrs <= -1] = -0.99999; diag(corrs) = 1
  omega = S %*% corrs %*% S
  
  fit.x = omega[1:Pz,Pw] %*% solve(omega[1:Pz,1:Pz]) %*% omega[1:Pz,Pw]
  fit.y = omega[1:Pz,P] %*% solve(omega[1:Pz,1:Pz]) %*% omega[1:Pz,P]  
  if (omega[Pw,Pw] <= fit.x) omega[Pw,Pw] = fit.x / 0.99999 #scale var up to have R2 slightly below 1
  fit.xy = omega[1:Pw,P] %*% solve(omega[1:Pw,1:Pw]) %*% omega[1:Pw,P]
  if (omega[P,P] <= fit.xy) omega[P,P] = fit.xy / 0.99999 #scale var up to have R2 slightly below 1
  
  p.cov = omega[Pw,P] - t(omega[1:Pz,Pw]) %*% solve(omega[1:Pz,1:Pz]) %*% omega[1:Pz,P]
  p.vars = c(omega[Pw,Pw] - fit.x, omega[P,P] - fit.y)
  out[1] = ifelse(all(p.vars > 0), p.cov/sqrt(prod(p.vars)), NA)
  
  draws = tryCatch(matrixsampling::rwishart(n.iter, K, Sigma=sigma/K, Theta=omega), error=function(e){return(NULL)}) 
  if (!is.null(draws)) {
    est = apply(draws, 3, estimate.pcor, sigma)  
    out[2:3] = quantile(est, c(0.025, 0.975), na.rm=T)    
    out[out < -1] = -1; out[out > 1] = 1
  }
  return(round(out,5))
}   



integral.p <- function (integral.func, K, omega, sigma, min.iter = 10000, adap.thresh = c(1e-04,  1e-06)) {
  tot.iter = min.iter * 10^(0:length(adap.thresh))
  adap.thresh = c(adap.thresh, 0)
  p = 1
  curr.iter = 0
  for (i in 1:length(tot.iter)) {
    add.iter = tot.iter[i] - curr.iter
    add.p = integral.func(K, omega, sigma, n.iter = add.iter)
    p = (curr.iter * p + add.iter * add.p)/tot.iter[i]
    curr.iter = tot.iter[i]
    if (all(is.na(p)) || all(p[!is.na(p)] >= adap.thresh[i])) 
      break
  }
  return(p)
}




pcov.integral <- function(K, omega, sigma, n.iter = 1000, add.reverse = T, xy.index = NULL, z.index = NULL) {
  P = dim(omega)[1]
  if (P <= 2) 
    return(NA)
  if (is.null(xy.index)) 
    xy.index = 1:2
  if (is.null(z.index)) 
    z.index = which(!(1:P %in% xy.index))
  index = check.index(xy.index, z.index, dim(omega)[1])
  if (is.null(index)) 
    return(NA)
  if (!add.reverse) {
    omega = omega[index, index]
    sigma = sigma[index, index]
    Pw = P - 1
    Pz = Pw - 1
    fit.x = omega[1:Pz, Pw] %*% solve(omega[1:Pz, 1:Pz]) %*% 
      omega[1:Pz, Pw]
    fit.y = omega[1:Pz, P] %*% solve(omega[1:Pz, 1:Pz]) %*% 
      omega[1:Pz, P]
    if (omega[Pw, Pw] <= fit.x) 
      omega[Pw, Pw] = fit.x/0.99999
    fit.xy = omega[1:Pw, P] %*% solve(omega[1:Pw, 1:Pw]) %*% 
      omega[1:Pw, P]
    if (omega[P, P] <= fit.xy) 
      omega[P, P] = fit.xy/0.99999
    var.y = as.numeric(sigma[P, P] - sigma[P, 1:Pw] %*% solve(sigma[1:Pw, 
                                                                    1:Pw]) %*% sigma[1:Pw, P])/K^2
    sig.xys = solve(sigma[1:Pw, 1:Pw]) %*% sigma[1:Pw, P]/K
    gamma.parts = list(x = solve(omega[1:Pz, 1:Pz]) %*% omega[1:Pz, 
                                                              Pw], y = solve(omega[1:Pz, 1:Pz]) %*% omega[1:Pz, 
                                                                                                          P], fit.x = fit.x * K, dw.dz.gamma = c(omega[1:Pz, 
                                                                                                                                                       P], omega[1:Pz, Pw] %*% solve(omega[1:Pz, 1:Pz]) %*% 
                                                                                                                                                   omega[1:Pz, P]) * K)
    sigma.use = matrix(0, P + Pw, P + Pw)
    sigma.use[1:Pw, 1:Pw] = sigma[1:Pw, 1:Pw]
    theta = matrix(0, P + Pw, P + Pw)
    theta[1:Pz + Pw, 1:Pz + Pw] = K * omega[1:Pz, 1:Pz]
    diag(theta)[Pw + Pw:P] = K * c(omega[Pw, Pw] - fit.x, 
                                   omega[P, P] - fit.y)
    draws = matrixsampling::rwishart(n.iter, K, Sigma = sigma.use, 
                                     Theta = theta)
    params = apply(draws, 3, pcov.cond.stats, K, sigma, gamma.parts, 
                   sig.xys, var.y)
    pcov.obs = omega[Pw, P] - t(omega[1:Pz, Pw]) %*% solve(omega[1:Pz, 
                                                                 1:Pz]) %*% omega[1:Pz, P]
    p = conditional.norm(pcov.obs, params[1, ], params[2, 
    ])
    return(p)
  }
  else {
    p1 = pcov.integral(K, omega, sigma, n.iter/2, add.reverse = F, 
                       xy.index = xy.index, z.index = z.index)
    p2 = pcov.integral(K, omega, sigma, n.iter/2, add.reverse = F, 
                       xy.index = rev(xy.index), z.index = z.index)
    return((p1 + p2)/2)
  }
}


filter.params <- function (data, locus.id, params, param.lim = 1.25, multreg = F) {
  out.of.bounds = (abs(data[params[1]]) > abs(param.lim)) | is.nan(as.numeric(data[params[1]]))
  if (any(out.of.bounds)) {
    if (!multreg) {
      phen.cols = c("phen1", "phen2")
    }
    else {
      phen.cols = c("predictors", "outcome")
    }
    message("Warning: Estimates too far out of bounds (+-", 
            param.lim, ") for phenotype(s) ", paste0(paste0(data[phen.cols[1]][out.of.bounds], 
                                                            " ~ ", data[phen.cols[2]][out.of.bounds], " (", 
                                                            signif(data[params[1]][out.of.bounds], 4), ")"), 
                                                     collapse = ", "), " in locus ", locus.id, ". Values will be set to NA. To change this threshold, modify the 'param.lim' argument")
    if (!multreg) {
      data[out.of.bounds, params] = NA
    }
    else {
      data[, params] = NA
    }
  }
  return(data)
}


cap <- function (values, lim = c(-1, 1)) {
  for (i in 1:length(values)) {
    if (is.na(values[i])) 
      next
    if (values[i] > max(lim)) {
      values[i] = max(lim)
    }
    else if (values[i] < min(lim)) {
      values[i] = min(lim)
    }
  }
  return(values)
}