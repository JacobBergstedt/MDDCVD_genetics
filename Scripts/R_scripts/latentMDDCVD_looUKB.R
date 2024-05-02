library(GenomicSEM)
library(tidyverse)
library(data.table)
library(lavaan)

# doesn't work with singularity, move there before submitting job
# setwd("/cluster/p/p33/cluster/users/jacobb/MDDCVD_sumstats")

args = commandArgs(trailingOnly=TRUE)


# Prepare the munged traits for CVD 
traits <- c("Data/Munged_sumstats/CAD_sumstats_munged.sumstats",
            "Data/Munged_sumstats/HF_sumstats_munged.sumstats",
            "Data/Munged_sumstats/PAD_sumstats_munged.sumstats",
            "Data/Munged_sumstats/Stroke_sumstats_munged.sumstats",
			      "Data/Munged_sumstats/MDD_PGC3_looUKB.sumstats.gz")
			
sample.prev <- c(0.215,0.0484,0.0269,0.129,0.21)
population.prev <- c(0.080,0.0484,0.074,0.076,0.16)
ld <- "Data/REF/eur_w_ld_chr"
wld <- "Data/REF/eur_w_ld_chr"
trait.names<-c("CAD","HF","PAD","Stroke","MDD_looUKB")

LDSCoutput <- ldsc(traits, sample.prev, population.prev, ld, wld, trait.names)


# Estimate model fit for hierarchical factor model

model <- 'CVD =~ 1*CAD + HF + PAD + Stroke
CVDMDD =~ 1*CVD + MDD_looUKB
CVD ~~ 0*CVD
'

EFA_CVD_MDD_looUKB <- usermodel(covstruc=LDSCoutput,
  estimation="DWLS",
  model=model)

EFA_CVD_MDD_looUKB

# Save the results

saveRDS(EFA_CVD_MDD_looUKB, "Results/gSEM_joelle/EFA_CVD_MDD_looUKB.rds")

# Prep MDD sumstats for processing

MDD <- fread("Data/Sumstats/daner_pgc_mdd_ex.23aME.loo.no.UKBB.gz")
MDD <- MDD %>% 
  rename(N_CAS=Nca, N_CON=Nco, RSID=SNP, POS=BP) %>% 
  mutate(B= log10(OR), N=N_CAS+N_CON, EAF=(FRQ_A_357636+FRQ_U_1281936)/N) %>% 
  select(CHR, RSID, POS, A1, A2, P, SE, N, B, EAF, N_CAS, N_CON, N)
  
fwrite(MDD, "Data/Sumstats/MDD_looUKB_sumstats", sep="\t")

# Make dataframe with sumstats for GWAS

files=c("Data/Sumstats/CAD.sumstats.gsem",
"Data/Sumstats/HF.sumstats.gsem",
"Data/Sumstats/PAD.sumstats.gsem",
"Data/Sumstats/STR.sumstats.gsem",
"Data/Sumstats/MDD_looUKB_sumstats")

ref= "Data/REF/reference.1000G.maf.0.005.txt"
trait.names=c("CAD","HF", "PAD", "Stroke", "MDD_looUKB")
se.logit=c(T,T,T,T,T)
maf.filter=0.01

p_sumstats <- sumstats(files=files,ref=ref,trait.names=trait.names,se.logit=se.logit, OLS=NULL,linprob=NULL,N=NULL,maf.filter=maf.filter, keep.indel=FALSE,parallel=F,cores=NULL)

saveRDS(p_sumstats, "Results/gSEM_joelle/EFA_CVD_loUKB_sumstats.rds")

p_sumstats <- p_sumstats %>% filter(CHR==args[1])

model <- 'CVD =~ 1*CAD + HF + PAD + Stroke
CVDMDD =~ 1*CVD + MDD_looUKB
CVD ~~ 0*CVD

CVDMDD ~ SNP
'

# Run the GWAS per chromosome

CVD_MDD_looUKB_GWAS <- userGWAS(covstruc = LDSCoutput, SNPs = p_sumstats, estimation = "DWLS", model = model, printwarn = TRUE, cores = 32, toler = FALSE, SNPSE = FALSE, parallel = TRUE)

saveRDS(CVD_MDD_looUKB_GWAS, paste0("Results/gSEM_joelle/CVD_MDD_looUKB_GWAS_",args[1],".rds"))

# Merge chromosome results 

files <- list.files("Results/gSEM_joelle/", pattern=glob2rx("CVD_MDD_looUKB_GWAS_*.rds"), full.names=T)

GWAS <- list()
GWAS2 <- list()

for (f in 1:length(files)){
GWAS[[f]] <- readRDS(files[f])
GWAS2[[f]] <- rbindlist(GWAS[[f]]) %>% 
  filter(lhs=="CVDMDD" & rhs=="SNP" & warning==0 & error==0)
full <- rbindlist(GWAS2)
fwrite(full, "Results/gSEM_joelle/CVD_MDD_looUKB_GWAS.txt", sep=" ", na=NA, quote=F)
}

# Estimate sample size

# restrict to MAF of 40% and 10%

GWAS3 <- subset(full, full$MAF <= .4 & full$MAF >= .1)

#calculate expected sample size (N_hat)

N_hat <- mean(1 / ((2 * GWAS3$MAF * (1 - GWAS3$MAF)) * GWAS3$SE^2))

write.csv(N_hat, "Results/gSEM_joelle/CVD_GWAS_looUKB_N.csv")
