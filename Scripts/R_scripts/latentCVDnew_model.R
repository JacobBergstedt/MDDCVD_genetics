library(GenomicSEM)
library(tidyverse)
library(data.table)
library(lavaan)
library(semPlot)
library(qqman)
library(ggplot2)

# setwd("/cluster/p/p33/cluster/users/jacobb/MDDCVD_sumstats")

# Now we use Als as MDD sumstats
# Jacob created ldscoutput with ldsc_gSEM.R

LDSCoutput <- readRDS("Data/Genetic_covariance_matrices/genetic_covariance_Als_2023.rds")

model_common <- '
CVD =~ 1*CAD_2022 + Stroke_Gigastroke_2022 + PAD_2021_V2 + HF_2019_V2
'


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
                                                         ">", "==", "|", "<", ">"), ]
       semModel@Vars <- data.frame(name = c(varNames, factNames),
                                   manifest = c(varNames, factNames) %in% varNames, exogenous = NA,
                                   stringsAsFactors = FALSE)
       semModel@ObsCovs <- list()
       semModel@ImpCovs <- list()
       semModel@Computed <- FALSE
       semModel@Original <- list(object)
       return(semModel)

}

forplot <- usermodel(covstruc=LDSCoutput,
                         estimation="DWLS",
                         model=model_common)


forplot$results <- forplot$results %>%
	mutate(lhs = ifelse(lhs == "Stroke_Gigastroke_2022", "STR", lhs)) %>%
	mutate(lhs = ifelse(lhs == "HF_2019_V2", "HF", lhs)) %>%
	mutate(rhs = ifelse(rhs == "Stroke_Gigastroke_2022", "STR", rhs)) %>%
	mutate(rhs = ifelse(rhs == "HF_2019_V2", "HF", rhs))

png("Results/gSEM_joelle/semplot_CVD.png", units="in", width=5, height=5, res=300)
semPaths(semPlotModel_GSEM(forplot),layout="tree2",whatLabels = "est")
dev.off()




