library(quantro)
library(quantroSim)
library(minfi)
library(genefilter)
library(preprocessCore)
library(doParallel)
library(abind)

source("quantro-functions.R")

###########################################
### Simulation study: DNA Methylation: ROC Plot
###########################################

# Set essential simulation scenario parameters
nMol = 1e6
nMethProbes = 1e4
nS = 5
nG = 2
totS = nS * nG
N = 1000
pDiffs = runif(N, 0.001, 0.50) 
nDiffs = nMethProbes * pDiffs
siga = sigb = 1 * diag(totS)
nCores = 10

acomb <- function(...) abind(..., along=3)
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
                                   pDiff = pDiffs[i], pUp = 0.80, verbose = FALSE)
    
    trueDiffProbes = array(0, dim = nMethProbes)
    trueDiffProbes[methTruth$methDiffInd] <- 1
    if(nG==2){ trueDiff = methTruth$methRange[,1] - methTruth$methRange[,2] }
    
    sim <- simulateMeth(methTruth, meth.platform = "methArrays", nSamps = nS,  
                        nMol = nMol, verbose = FALSE, siga = siga, sigb = sigb)
    
    Mset1 <- preprocessMeth(simulateMethObject = sim, method = "quantile")
    dmpPicksMset1 <- rowttests(as.matrix(Mset1$beta), sim$pd$Group)
    pvalMset1 <- dmpPicksMset1$p.value
    
    Mset2 = preprocessMeth(simulateMethObject = sim, method = "quantro", 
                           B = 100, alpha = 0.05) 
    dmpPicksMset2 <- rowttests(as.matrix(Mset2$beta), sim$pd$Group)
    pvalMset2 <- dmpPicksMset2$p.value
    
    pBH1 <- p.adjust(pvalMset1, method = "BH")
    pBH2 <- p.adjust(pvalMset2, method = "BH")
    
    x = data.frame("truth" = trueDiffProbes, "pBH1" = pBH1, "pBH2" = pBH2)
    o.pBH1 <- order(x$pBH1)
    o.pBH2 <- order(x$pBH2)
    
    z = data.frame(array(0, dim = c(nMethProbes, 8)))
    colnames(z) = c("TP1", "TP2", "FP1", "FP2", "FN1", "FN2", "TN1", "TN2")
    for(j in seq_len(nMethProbes)){
        z[j,1] = length(which(x[head(o.pBH1, j),]$truth == 1))  
        z[j,2] = length(which(x[head(o.pBH2, j),]$truth == 1))   
        z[j,3] = length(which(x[head(o.pBH1, j),]$truth == 0))  
        z[j,4] = length(which(x[head(o.pBH2, j),]$truth == 0)) 
        z[j,5] = length(which(x[tail(o.pBH1, -j),]$truth == 1))  
        z[j,6] = length(which(x[tail(o.pBH2, -j),]$truth == 1))    
        z[j,7] = length(which(x[tail(o.pBH1, -j),]$truth == 0))  
        z[j,8] = length(which(x[tail(o.pBH2, -j),]$truth == 0))
    }
    
    TPR1 = z$TP1 / (z$TP1 + z$FN1)
    TPR2 = z$TP2 / (z$TP2 + z$FN2)
    FPR1 = z$FP1 / (z$FP1 + z$TN1)
    FPR2 = z$FP2 / (z$FP2 + z$TN2)
 
    return(list(TPR1, TPR2, FPR1, FPR2))
}

TPR1 <- simResults[[1]]
TPR2 <- simResults[[2]]
FPR1 <- simResults[[3]]
FPR2 <- simResults[[4]]

dat = data.frame("TPR" = c(rowMeans(TPR1), rowMeans(TPR2)), 
                 "FPR" = c(rowMeans(FPR1), rowMeans(FPR2)), 
                 "Normalization" = c(rep("quantile", nMethProbes), 
                                     rep("quantro", nMethProbes)))

# partial AUC function
pAUC <- function(fpr, tpr, lower=0, upper=1){
    ininterval <- which( fpr>=lower & fpr<=upper )
    n <- length(ininterval)
    
    delta <- diff(fpr[ininterval])
    av <- (tpr[ininterval][-1]+tpr[ininterval][-n])/2
    sum(delta*av)
}

x = dat %>% group_by(Normalization) %>% summarise(pAUC = pAUC(FPR, TPR, upper = 0.25))
x$pAUC <- round(x$pAUC, 3)

# Plot ROC curves
ggplot(dat, aes(x=FPR, y=TPR, color = Normalization)) + 
    geom_line(size =1) + xlab("False Positive Rate") + 
    ylab("True Positive Rate") + xlim(0, 0.25) +
    labs(title="Proportion of Differences Between Groups\n(pDiff ~ Uniform[0, 0.50])") + 
    scale_colour_discrete(breaks=c("quantile", "quantro"),
                          labels=c(paste0("quantile (pAUC: ", x$pAUC[1],")"), 
                                   paste0("quantro (pAUC: ", x$pAUC[2],")"))) +
    theme(axis.text=element_text(size=16, color = "black"),
          axis.title=element_text(size=16,face="bold"), 
          title = element_text(size = 16), 
          legend.title = element_text(size = 20),
          legend.text = element_text(size = 20),
          legend.justification=c(1,0), legend.position=c(1,0))



