library(quantro)
library(quantroSim)
library(minfi)
library(genefilter)
library(preprocessCore)
library(doParallel)
library(abind)
library(ggplot2)
library(dplyr)

source("quantro-functions.R")

###########################################
### Simulation study: FDR plots
###########################################

# Set essential simulation scenario parameters
nMol = 1e6
nMethProbes = 450*1e3 # simulate a 450k array
nS = 5
nG = 2
totS = nS * nG
pDiffs = 0.25
nDiffs = nMethProbes * pDiffs
siga = sigb = 1 * diag(totS)
N = 100
nCores = 5

registerDoParallel(cores = nCores)
workers <- getDoParWorkers()
backend <- getDoParName()
version <- getDoParVersion()

comb <- function(x, a) {
    lapply(seq_along(a), function(ii){
        abind::abind(x[[ii]], a[[ii]], along = 2) })
}

simResults <- foreach(i = 1:N, .combine='comb') %dopar% { 
    methTruth <- simulateMethTruth(nProbes = nMethProbes, nGroups = nG, 
                                   pDiff = pDiffs, pUp = 0.80, verbose = FALSE)
    
    trueDiffProbes = array(0, dim = nMethProbes)
    trueDiffProbes[methTruth$methDiffInd] <- 1
    
    sim <- simulateMeth(methTruth, meth.platform = "methArrays", nSamps = nS,  
                        nMol = nMol, verbose = FALSE, siga = siga, sigb = sigb)
    
    Mset1 <- preprocessMeth(simulateMethObject = sim, method = "quantile")
    pvalMset1 <- rowttests(as.matrix(Mset1$beta), sim$pd$Group)$p.value
    
    Mset2 = preprocessMeth(simulateMethObject = sim, method = "quantro", 
                           B = 100, alpha = 0.05) 
    pvalMset2 <- rowttests(as.matrix(Mset2$beta), sim$pd$Group)$p.value
    
    pBH1 <- p.adjust(pvalMset1, method = "BH")
    pBH2 <- p.adjust(pvalMset2, method = "BH")
    
    x = data.frame("truth" = trueDiffProbes, "pBH1" = pBH1, "pBH2" = pBH2)
    o.pBH1 <- order(x$pBH1)
    o.pBH2 <- order(x$pBH2)
    
    z1 = z2 = array(0, dim = nDiffs)
    for(j in seq_len(nDiffs)){
        z1[j] = length(which(x[head(o.pBH1, j),]$truth == 0))  
        z2[j] = length(which(x[head(o.pBH2, j),]$truth == 0))  
    }
    return(list(z1, z2))      
}

z1 <- simResults[[1]]
z2 <- simResults[[2]]

dat = data.frame("mean" = c(rowMeans(z1), rowMeans(z2)), 
                 "se" = c(apply(z1, 1, sd) / sqrt(N), apply(z2, 1, sd) / sqrt(N)),
                 "Normalization" = c(rep("quantile", nDiffs), rep("quantro", nDiffs)),
                 "nSelected" = rep(seq_len(nDiffs),2))         
x = dat %>% group_by(Normalization) %>% summarise(max = round(mean[max(nSelected)]))

# Create plot using ggplot2
ggplot(dat, aes(x=nSelected, y=mean, color = Normalization)) + 
    geom_line(size = 1) + 
    xlab("Number of Top Probes Selected") + 
    ylab("Number of False Discoveries") + 
    labs(title=paste0("nProbes=", nMethProbes, ", pDiff=", pDiffs)) + 
    scale_colour_discrete(breaks=c("quantile", "quantro"),
                          labels=c(paste0("quantile (# false discoveries: ", x$max[1],")"), 
                                   paste0("quantro (# false discoveries: ", x$max[2],")"))) +
    theme(axis.text=element_text(size=16, color = "black"),
          axis.title=element_text(size=16,face="bold"), 
          title = element_text(size = 16), 
          legend.title = element_text(size = 20),
          legend.text = element_text(size = 20),
          legend.justification=c(0,1), legend.position=c(0,1))



