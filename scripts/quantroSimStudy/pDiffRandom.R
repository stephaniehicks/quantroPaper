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
### Simulation study: Bias and MSE
###########################################

# Set essential simulation scenario parameters
nMol = 1e6
nMethProbes = 1e4
nS = 5
nG = 2
totS = nS * nG
N = 1000
nCores = 5
pDiffs = runif(N, 0.001, 0.08)
nDiffs = nMethProbes * pDiffs
pUp = 0.80
magTechVar <- c(.1, 1)
magTechVarOpt = c(1, 2)
quantroAlpha = c(seq(0.01, 0.04, by = 0.01), seq(0.05, 0.95, by = 0.05))
lenP <- length(quantroAlpha)
nMethods <- 2 + lenP

#### Low technical variation
simResults <- compareMethodsRandom(nMethProbesObject = nMethProbes, 
                    pDiffObject = pDiffs, pUp = pUp, 
                    meth.platform = "methArrays", 
                    nGroups = nG, nSamps = nS, nMol = nMol, nCores = nCores, 
                    B = NULL, quantroAlpha = quantroAlpha,
                    siga = magTechVar[1]*diag(totS), 
                    sigb = magTechVar[1]*diag(totS), 
                    sigOpt = magTechVarOpt[1]*diag(totS))

simDiffs <- abind(simResults, along = 3)
outMSE <- calculateMSERandom(simDiffs, nGroups = nG)
MSEMat <- abind(outMSE, along = 2) 
colnames(MSEMat) <- c("Mean", "Bias", "Bias2", "Var", "MSE")

outLoss <- calculateLossRandom(simDiffs, nGroups = 2, topK = 10)
LossMat <- colMeans(outLoss$simLossMat)
OverlapMat <- colMeans(outLoss$overlapMat)


#### High technical variation
simResults <- compareMethodsRandom(nMethProbesObject = nMethProbes, 
                                   pDiffObject = pDiffs, pUp = pUp, 
                                   meth.platform = "methArrays", 
                                   nGroups = nG, nSamps = nS, nMol = nMol, nCores = nCores, 
                                   B = NULL, quantroAlpha = quantroAlpha,
                                   siga = magTechVar[2]*diag(totS), 
                                   sigb = magTechVar[2]*diag(totS), 
                                   sigOpt = magTechVarOpt[2]*diag(totS))

simDiffs <- abind(simResults, along = 3)
outMSE <- calculateMSERandom(simDiffs, nGroups = nG)
MSEMat <- abind(outMSE, along = 2) 
colnames(MSEMat) <- c("Mean", "Bias", "Bias2", "Var", "MSE")

outLoss <- calculateLossRandom(simDiffs, nGroups = 2, topK = 10)
LossMat <- colMeans(outLoss$simLossMat)
OverlapMat <- colMeans(outLoss$overlapMat)


