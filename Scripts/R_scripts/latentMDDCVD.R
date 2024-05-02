library(GenomicSEM)
library(tidyverse)
library(data.table)
# library(lavaan)
# library(semPlot)
# library(qqman)
# library(ggplot2)

# doesn't work with singularity, move there before submitting job
# setwd("/cluster/p/p33/cluster/users/jacobb/MDDCVD_sumstats")

# In this script we use Howard 2018 as measure for MDD
# Ncases=40675, Ncontrols=64643, prev=0.39

args = commandArgs(trailingOnly=TRUE)

##################################################################################################
# ###A common factor model underlying CVD and MDD
 

# Prepare the munged traits for CVD 
traits <- c("Data/Munged_sumstats/CAD_sumstats_munged.sumstats",
"Data/Munged_sumstats/HF_sumstats_munged.sumstats",
"Data/Munged_sumstats/PAD_sumstats_munged.sumstats",
"Data/Munged_sumstats/Stroke_sumstats_munged.sumstats",
"Data/Munged_sumstats/MDD_2018_Howard_sumstats_munged.sumstats")
			
sample.prev <- c(0.215,0.0484,0.0269,0.129,0.111,0.39)
population.prev <- c(0.080,0.0484,0.074,0.076,0.111,0.16)
ld <- "/REF/ldsc/eur_w_ld_chr"
wld <- "/REF/ldsc/eur_w_ld_chr"
trait.names<-c("CAD","HF","PAD","STR","MDD2018")

LDSCoutput <- ldsc(traits, sample.prev, population.prev, ld, wld, trait.names)

LDSCcov <- LDSCoutput$S
k <- nrow(LDSCoutput$S)
LDSCcov_SE <- matrix(0, k, k)
LDSCcov_SE[lower.tri(LDSCcov_SE,diag=TRUE)] <-sqrt(diag(LDSCoutput$V))


## Estimate model fit for common factor model

EFA_CVD_MDD <- commonfactor(covstruc = LDSCoutput, estimation="DWLS")
EFA_CVD_MDD

## Save the results

saveRDS(LDSCoutput, "Results/gSEM_joelle/EFA_CVD_MDD2018_ldscoutput.rds")
saveRDS(LDSCcov_SE, "Results/gSEM_joelle/EFA_CVD_MDD2018_ldscoutput_SE.rds")
saveRDS(EFA_CVD_MDD_1, "Results/gSEM_joelle/EFA_CVD_MDD2018.rds")

## Create a plot
## https://groups.google.com/g/genomic-sem-users/c/h250OMkp_CM
## Function from Nicola Pirastu
semPlotModel_GSEM=function(gsem.object=GWISoutput , est.label="STD_All"){
       object=gsem.object$results
       object$free=0
       numb=1:length(which(object$op!="~~"))
       object$free[which(object$op!="~~")]=numb
       varNames <- lavaanNames(object, type = "ov")
       factNames <- lavaanNames(object, type = "lv")
       factNames <- factNames[!factNames %in% varNames]
       n <- length(varNames)
       k <- length(factNames)
       if (is.null(object$label))
       object$label <- rep("", nrow(object))
       semModel <- new("semPlotModel")
       object$est <- object[,est.label]
       if (is.null(object$group))
       object$group <- ""
       semModel@Pars <- data.frame(label = object$label, lhs = ifelse(object$op ==
       "~" | object$op == "~1", object$rhs, object$lhs), edge = "--",
       rhs = ifelse(object$op == "~" | object$op == "~1", object$lhs,
       object$rhs), est = object$est, std = NA, group = object$group,
       fixed = object$free==0, par = object$free, stringsAsFactors = FALSE)
       semModel@Pars$edge[object$op == "~~"] <- "<->"
       semModel@Pars$edge[object$op == "~*~"] <- "<->"
       semModel@Pars$edge[object$op == "~"] <- "~>"
       semModel@Pars$edge[object$op == "=~"] <- "->"
       semModel@Pars$edge[object$op == "~1"] <- "int"
       semModel@Pars$edge[grepl("\\|", object$op)] <- "|"
       semModel@Thresholds <- semModel@Pars[grepl("\\|", semModel@Pars$edge),
       -(3:4)]
       semModel@Pars <- semModel@Pars[!object$op %in% c(":=", "<",
                                                        # ">", "==", "|", "<", ">"), ]
       semModel@Vars <- data.frame(name = c(varNames, factNames),
                                   manifest = c(varNames, factNames) %in% varNames, exogenous = NA,
                                   stringsAsFactors = FALSE)
       semModel@ObsCovs <- list()
       semModel@ImpCovs <- list()
       semModel@Computed <- FALSE
       semModel@Original <- list(object)
       return(semModel)

}

model <- 'CVD =~ 1*MDD + HF + PAD + STR + CAD'

forplot <- usermodel(covstruc=LDSCoutput,
                         estimation="DWLS",
                         model=model)


png("Results/gSEM_joelle/CVD_MDD.png", units="in", width=5, height=5, res=300)
semPaths(semPlotModel_GSEM(forplot),layout="tree2",whatLabels = "est")
dev.off()




##################################################################################################
###A hierarchical factor model underlying CVD and MDD
 
LDSCoutput <- readRDS("Results/gSEM_joelle/EFA_CVD_MDD2018_ldscoutput.rds")

# Estimate model fit for hierarchical factor model
# Note we do not include AF due to it's low rg/ factor loading with the other CVDs

model <- 'CVD =~ 1*CAD + HF + PAD + STR
CVDMDD =~ 1*MDD2018 + CVD
CVD ~~ 0*CVD
'

EFA_CVD_MDD_2 <- usermodel(covstruc=LDSCoutput, estimation="DWLS",  model=model)
EFA_CVD_MDD_2

## Save the results

saveRDS(EFA_CVD_MDD_2, "Results/gSEM_joelle/EFA_CVD_MDD2018_2.rds")

png("Results/gSEM_joelle/CVD_MDD2018.png", units="in", width=5, height=5, res=300)
semPaths(semPlotModel_GSEM(EFA_CVD_MDD_2),layout="tree2",whatLabels = "est")
dev.off()















##################################################################################################
### Prep the clean sumstats files for the GWASs
### The software doesn't like having 2 effect columns, so remove Z-column

CAD <- fread("sumstats/CAD.sumstats")
CAD <- CAD %>%
select(-Z)
fwrite(CAD, "sumstats/CAD.sumstats.gsem", sep="\t")

HF <- fread("sumstats/HF.sumstats")
HF <- HF %>%
select(-Z)
fwrite(HF, "sumstats/HF.sumstats.gsem", sep="\t")

PAD <- fread("sumstats/sumstats_nadine_cleaned/PAD_cleaned_GRCh37.sumstats")
PAD <- PAD %>%
select(-Z) %>%
rename(c(A1=EffectAllele, A2=OtherAllele))
fwrite(PAD, "sumstats/PAD.sumstats.gsem", sep="\t")

STR <- fread("sumstats/sumstats_nadine_cleaned/stroke_cleaned_GRCh37.sumstats") %>%
select(-Z) %>%
rename(c(A1=EffectAllele, A2=OtherAllele))
fwrite(STR, "sumstats/STR.sumstats.gsem", sep="\t")

MDD <- fread("sumstats/sumstats_nadine_cleaned/MDD_cleaned_GRCh37.sumstats") %>%
select(-Z) %>%
rename(c(A1=EffectAllele, A2=OtherAllele))
fwrite(MDD, "sumstats/MDD.sumstats.gsem", sep="\t")

IL6 <- fread("sumstats/IL6.sumstats")
IL6 <- IL6 %>%
select(-Z)
fwrite(IL6, "sumstats/IL6.sumstats.gsem", sep="\t")

MDD2018 <- fread("Data/Sumstats/MDD_2018_Howard_sumstats")
MDD2018 <- MDD2018 %>%
select(-Z) %>%
filter(as.numeric(CHR)<23)
fwrite(MDD2018, "Data/Sumstats/MDD_2018_Howard_sumstats.gsem", sep="\t")

##################################################################################################
###GWAS on common factor model underlying CVD 

LDSCoutput <- readRDS("genomicSEMresults/EFA_CVD_1_ldscoutput.rds")

files=c("sumstats/CAD.sumstats.gsem",
"sumstats/HF.sumstats.gsem",
"sumstats/PAD.sumstats.gsem",
"sumstats/STR.sumstats.gsem")
			
ref= "reference.1000G.maf.0.005.txt"
trait.names=c("CAD","HF", "PAD", "STR")
se.logit=c(T,T,T,T)
maf.filter=0.01

p_sumstats <- sumstats(files=files,ref=ref,trait.names=trait.names,se.logit=se.logit,
                      # OLS=NULL,linprob=NULL,prop=NULL,N=NULL,maf.filter=maf.filter,
                      # keep.indel=FALSE,parallel=F,cores=NULL)
					  
saveRDS(p_sumstats, "genomicSEMresults/EFA_CVD_sumstats.rds")

p_sumstats <- readRDS("genomicSEMresults/EFA_CVD_sumstats.rds")

CVD_GWAS <-commonfactorGWAS(covstruc = LDSCoutput, SNPs = p_sumstats, estimation = "DWLS", cores=32, toler = FALSE, SNPSE = FALSE, parallel = T, GC="standard")

saveRDS(CVD_GWAS, "genomicSEMresults/CVD_GWAS.rds")

# ### Estimate sample size

# #restrict to MAF of 40% and 10%
GWAS2<-subset(CVD_GWAS, CVD_GWAS$MAF <= .4 & CVD_GWAS$MAF >= .1)

# #calculate expected sample size (N_hat)
N_hat<-mean(1/((2*GWAS2$MAF*(1-GWAS2$MAF))*GWAS2$se_c^2))
N_hat

write.csv(N_hat, "CVD_GWAS_N.csv")

# #export sumstats for SNP effects on the common latent factor

GWAS3 <- CVD_GWAS %>%
filter(lhs=="F1" & op=="~" & rhs=="SNP" & warning==0 & fail==0)

GWAS4 <- GWAS3 %>%
filter(Q_pval>0.05)

fwrite(GWAS3, "CVD_GWAS.txt")
fwrite(GWAS4, "CVD_GWAS_Q0.05.txt")


	

##################################################################################################
###GWAS on common & bifactor factor model underlying CVD and MDD

LDSCoutput <- readRDS("Results/gSEM_joelle/EFA_CVD_MDD2018_ldscoutput.rds")

files=c("Data/Sumstats/CAD.sumstats.gsem",
"Data/Sumstats/HF.sumstats.gsem",
"Data/Sumstats/PAD.sumstats.gsem",
"Data/Sumstats/STR.sumstats.gsem",
"Data/Sumstats/MDD_2018_Howard_sumstats.gsem"
)
			
ref= "Data/Sumstats/reference.1000G.maf.0.005.txt"
trait.names=c("CAD","HF", "PAD", "STR", "MDD2018")
se.logit=c(T,T,T,T,T)
maf.filter=0.01

p_sumstats <-sumstats(files=files,ref=ref,trait.names=trait.names,se.logit=se.logit,
OLS=NULL,linprob=NULL,prop=NULL,N=NULL,maf.filter=maf.filter,
keep.indel=FALSE,parallel=F,cores=NULL)

saveRDS(p_sumstats, "Results/gSEM_joelle/EFA_CVD_MDD2018_sumstats.rds")

p_sumstats <- readRDS("genomicSEMresults/EFA_CVD_MDD_sumstats.rds")

CVD_MDD_GWAS <-commonfactorGWAS(covstruc = LDSCoutput, SNPs = p_sumstats, estimation = "DWLS", cores=16, toler = FALSE, SNPSE = FALSE, parallel = T, GC="standard",MPI=FALSE)

saveRDS(CVD_MDD_GWAS, "genomicSEMresults/CVD_MDD_GWAS.rds")

GWAS <- readRDS("genomicSEMresults/CVD_MDD_GWAS.rds")

# ### Estimate sample size

# #restrict to MAF of 40% and 10%
GWAS2<-subset(GWAS, GWAS$MAF <= .4 & GWAS$MAF >= .1)

# #calculate expected sample size (N_hat)
N_hat<-mean(1/((2*GWAS2$MAF*(1-GWAS2$MAF))*GWAS2$SE^2))
N_hat

write.csv(N_hat, "CVD_MDD_GWAS_N.csv")

# #export sumstats for SNP effects on the common latent factor

GWAS3 <- GWAS %>%
	filter(lhs=="F1" & op=="~" & rhs=="SNP" & warning==0 & fail==0)
	
fwrite(GWAS3, "CVD_MDD_GWAS.txt")

# # #filter out SNPs with high Q

GWAS_Q <- GWAS3 %>%
filter(Q_pval>0.05)
	
# #Extract SNPs with some effect on MDD	
MDD <- fread("sumstats/MDD.sumstats.gsem")
MDDhits <- MDD %>%
rename(SNP=RSID) %>%
filter(P<0.05) %>%
mutate(CHR=as.numeric(CHR))
	
	
GWAS_Q2 <- GWAS_Q %>%
inner_join(MDDhits)
	
# #export homogeneous sumstats
fwrite(GWAS_Q2, "CVD_MDD_GWAS_filteredQ.txt")


##################################################################################################
###GWAS on bifactor factor model underlying CVD and MDD

# Full ldscoutput dataframe to experiment with: /cluster/p/p33/cluster/users/jacobb/MDDCVD_sumstats/Data/Genetic_covariance_matrices/genetic_covariance.rds
LDSCoutput <- readRDS("/cluster/p/p33/cluster/users/jacobb/MDDCVD_sumstats/Data/Genetic_covariance_matrices/genetic_covariance.rds")

LDSCoutput <- readRDS("Results/gSEM_joelle/EFA_CVD_MDD2018_ldscoutput.rds")
p_sumstats <- readRDS("Results/gSEM_joelle/EFA_CVD_MDD2018_sumstats.rds")

p_sumstats <- p_sumstats %>% filter(CHR==args[1])

model <- 'CVD =~ 1*CAD + HF + PAD + STR
CVDMDD =~ 1*CVD + MDD2018
CVD ~~ 0*CVD

CVDMDD ~ SNP
'

CVD_MDD2018_GWAS2 <- userGWAS(covstruc = LDSCoutput, SNPs = p_sumstats, estimation = "DWLS", model = model, printwarn = TRUE, cores = 32, toler = FALSE, SNPSE = FALSE, parallel = TRUE)

saveRDS(CVD_MDD2018_GWAS2, paste0("Results/gSEM_joelle/CVD_MDD_GWAS2_",args[1],".rds"))


##################################################################################################
###Heterogeneity for bifactor factor model underlying CVD and MDD

# First, re-run the GWAS with independent pathways instead of single common pathway

model <- 'CVD =~ 1*CAD + HF + PAD + STR
CVDMDD =~ 1*CVD + MDD2018
CVD ~~ 0*CVD

CAD + HF + PAD + STR + MDD2018 ~ SNP
'
	
CVD_MDD2018_GWAS_indep <- userGWAS(covstruc = LDSCoutput, SNPs = p_sumstats, estimation = "DWLS", model = model, printwarn = TRUE, cores = 32, toler = FALSE, SNPSE = FALSE, parallel = TRUE)

saveRDS(CVD_MDD2018_GWAS_indep, paste0("Results/gSEM_joelle/CVD_MDD2018_GWAS_indep_",args[1],".rds"))

# Merge chromosome results for independent pathway model

files <- list.files("genomicSEMresults", pattern=glob2rx("Results/gSEM_joelle/CVD_MDD2018_GWAS_indep_*rds"), full.names=T)

GWAS_indep <- list()
GWAS2_indep <- list()
for (f in 1:length(files)){
 GWAS_indep[[f]] <- readRDS(files[f])
 GWAS2_indep[[f]] <- rbindlist(GWAS_indep[[f]])
 full_indep <- rbindlist(GWAS2_indep)
 fwrite(full_indep, "Results/gSEM_joelle/CVD_MDD_GWAS_indep.txt", sep=" ")
}

# Merge chromosome results for common pathway model 

files <- list.files("Results/gSEM_joelle", pattern="CVD_MDD_GWAS2_", full.names=T)

GWAS <- list()
GWAS2 <- list()
for (f in 1:length(files)){
GWAS[[f]] <- readRDS(files[f])
GWAS2[[f]] <- rbindlist(GWAS[[f]]) %>% filter(lhs=="CVDMDD" & rhs=="SNP" & warning==0 & error==0)
full <- rbindlist(GWAS2)
fwrite(full, "Results/gSEM_joelle/CVD_MDD_GWAS_2.txt", sep=" ", na="NA")
}




# Calculate chisq difference per SNP (chisq is same for all paths in the model)

indep <- full_indep %>%
select(SNP, chisq, chisq_df) %>%
rename(chisq_indep=chisq) %>%
rename(chisq_df_indep=chisq_df) %>%
distinct()

GWAS2 <- full %>%
inner_join(indep) %>%
mutate(Q_chisq=chisq-chisq_indep) %>%
mutate(df=chisq_df-chisq_df_indep) %>%
mutate(Q_chisq_pval = pchisq(Q_chisq, df, lower.tail=FALSE))

GWAS3 <- GWAS2 %>%
filter(Q_chisq_pval>0.05)
	
Qsnps <- GWAS2 %>%
filter(Q_chisq_pval<5e-08 & Pval_Estimate<5e-08) %>%
select(SNP)

Qsnps <- as.character(Qsnps$SNP)
	
fwrite(GWAS2, "genomicSEMresults/CVD_MDD_GWAS_2.txt", sep=" ", na="NA")
fwrite(GWAS3, "genomicSEMresults/CVD_MDD_GWAS_2_Q0.05.txt", sep=" ", na="NA")

## Format results for single cell follow-up analyses (Arvid)

#CHR BP SNP A1 A2 FREQ B Z SE P N NCAS NCON INFO

formatted <- GWAS3 %>%
	select(CHR, BP, SNP, A1, A2, MAF, est, Z_Estimate, SE, Pval_Estimate) %>%
	mutate(N=333512)

names(formatted) <- c("CHR", "BP", "SNP", "A1", "A2", "FREQ", "B", "Z", "SE", "P", "N")

fwrite(formatted, "genomicSEMresults/CVD_MDD_GWAS_2_Q0.05_formatted.txt")

##################################################################################################

# Manhattan plot with highlight for heterogenous SNPs


GWAS2 <- fread("genomicSEMresults/CVD_MDD_GWAS_2.txt")

Qsnps <- GWAS2 %>%
filter(Q_chisq_pval<0.05 & Pval_Estimate<5e-08) %>%
select(SNP)

Qsnps <- as.character(Qsnps$SNP)

tiff("CVD_MDD_GWAS_manhattan.tiff", width=25, height=20, units="cm", res=300)
manhattan(GWAS2,
main = "Latent CVD-MDD",
chr="CHR",
bp="BP",
snp="SNP",
p="Pval_Estimate",
ylim=c(0,45),
cex = 1,
col = c("antiquewhite4", "deepskyblue4"),
highlight = Qsnps)
dev.off()

# Manhattan plot with CVD in the mirror

GWAS3 <- fread("CVD_MDD_GWAS_2_Q0.05.txt")
CVD <- fread("CVD_GWAS_Q0.05.txt")


tiff("CVD_MDD_mirror_manhattan.tiff",  width=25, height=30, units="cm", res=300)
par(mfrow=c(2,1))
par(mar=c(0,5,3,3))
manhattan(GWAS3,
p="Pval_Estimate",
ylim=c(0,35),
cex=1,
cex.lab=2.5,
font.lab=2,
font.axis=2,
cex.axis=1.6,
las=2,
col=c("antiquewhite4", "deepskyblue4"))
par(mar=c(5,5,3,3))
manhattan(CVD,
p="Pval_Estimate",
ylim=c(35,0),
cex=1.5,
cex.lab=2.5,
font.lab=2,
font.axis=2,
cex.axis=1.6,
las=2,
col=c("antiquewhite4", "deepskyblue4"),
xlab="",
xaxt="n")
dev.off()


