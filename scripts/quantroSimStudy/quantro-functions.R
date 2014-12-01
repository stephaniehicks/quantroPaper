#' @title Split it
#'
#' @description Split it
#' 
#' 
#' @param x
#'
#' @export 
splitit <- function(x){ split(seq(along=x),x) }


#' @title Preprocess the simulated DNA methylation data
#'
#' @description After simulating the DNA methylation data, preprocess the data using 
#' either no normalization, quantile normalization or normalization determined by quantro() in quantro package
#' 
#' @param simulateMethTruthObject Must be an object created from \code{simulateMethTruth}. 
#' @param method method to normalize arrays. Options include c("none", "quantile", "quantro"). Default is "none". 
#' @param B an optional argument. If \code{B} is not equal to NULL, 
#' then \code{B} represents the number of bootstrapped samples for 
#' permutation testing to assess statistical significance. 
#' @param alpha 
#' @param verbose
#'
#' @export
#' 
preprocessMeth <- function(simulateMethObject, method = "none", B = NULL, 
                           alpha = NULL, verbose = FALSE)
{
    if(exists("objectType", where = simulateMethObject)){
        if( !(simulateMethObject$typePlatform %in% c("methArrays")) ){
            stop("Platform currently not supported. Must be a platform list.meth.platforms().")
        }
        typePlatform <- simulateMethObject$typePlatform
        nProbes <- simulateMethObject$nProbes
        nSamps <- simulateMethObject$nSamps
        nGroups <- simulateMethObject$nGroups
        pd <- simulateMethObject$pd
        gID <- simulateMethObject$params$gID
        methDat <- simulateMethObject$meth
        unmethDat <- simulateMethObject$unmeth
        betaDat <- getBeta(getMethylSet(simulateMethObject), offset = 100)
        
    } else { stop("Must provide a DNA methylation object from simulateMeth().") }
    
    if(!(method %in% c("none", "quantile", "quantro"))){
        stop("Preprocessing method must be: none, quantile or quantro.")
    }
    
    if(method == "none"){
        if(verbose){ message("Preprocessing: no normalization.")}
        MsetObject <- list("objectType" = "simulateMethObject", "typePlatform" = typePlatform, "preprocessMethod" = method,
                           "nProbes" = nProbes, "nSamps" = nSamps, "nGroups" = nGroups, 
                           "pd" = pd, "meth" = methDat, "unmeth" = unmethDat, "beta" = betaDat)
    }
    
    if(method == "quantile"){
        if(verbose){ message("Preprocessing: quantile normalization.")}
        betaNorm = preprocessCore::normalize.quantiles(betaDat)
        colnames(betaNorm) <- pd$Sample_Name
        MsetObject <- list("objectType" = "simulateMethObject", "typePlatform" = typePlatform, "preprocessMethod" = method, 
                           "nProbes" = nProbes, "nSamps" = nSamps, "nGroups" = nGroups, 
                           "pd" = pd, "meth" = methDat, "unmeth" = unmethDat, 
                           "beta" = betaNorm)
    }
    
    if(method == "quantro"){
        if(verbose){ message("Preprocessing: quantro.")}
        betaFit <- quantro(betaDat, groupFactor = pd$Group, 
                           B = B, verbose = FALSE)
        
        if(quantroPvalPerm(betaFit) >= alpha){
            betaNorm <- normalize.quantiles(betaDat)
            colnames(betaNorm) <- pd$Sample_Name
            preproc = "quantile-quantro"
        } else { 
            betaNorm <- betaDat
            preproc = "none-quantro" 
        }
    
        MsetObject <- list("objectType" = "simulateMethObject", 
                           "typePlatform" = typePlatform, 
                           "preprocessMethod" = preproc, "nProbes" = nProbes, 
                           "nSamps" = nSamps, "nGroups" = nGroups, "pd" = pd, 
                           "quantroBeta" = betaFit, "meth" = methDat, 
                           "unmeth" = unmethDat, "beta" = betaNorm)
    }    
    return(MsetObject)
}


#' @title Get differentially methylated positions (DMPs) using a t-test
#'
#' @description Get differentially methylated positions (DMPs) using a t-test
#' 
#' @param simulateMethObject methArray object from simulateMethArrays() or preprocessMeth()
#' @param tTestAlpha
#' @param offset 
#'
#' @export
#' 
getDMPs <- function(simulateMethObject, tTestAlpha = 0.05, offset = 100)
{
    if(exists("objectType", where = simulateMethObject)){
        if( !(simulateMethObject$typePlatform %in% c("methArrays")) ){
            stop("Platform currently not supported. Must be a platform list.meth.platforms().")
        }
        nProbes = simulateMethObject$nProbes
        nGroups <- simulateMethObject$nGroups		
        pd <- simulateMethObject$pd
        betaDat <- getBeta(getMethylSet(simulateMethObject), offset = offset)
    } else { 
        stop("Must provide a DNA methylation object from simulateMeth() or preprocessMeth().") 
    }
    
    # t-tests
    tIndexes= splitit(pd$Group)
    tstatList = lapply(tIndexes, function(i) {
        x = rep(0,ncol(betaDat))
        x[i] = 1
        return(rowttests(betaDat, factor(x)))
    })
    
    if(nGroups == 2){
        tstatList <- tstatList[[2]]
        # pvalList <- tstatList[[2]][,"p.value"]
    }
    
    return(tstatList)
}	



compareMethodsRandom <- function(nMethProbesObject, pDiffObject, pUp = pUp,
                                 meth.platform = "methArrays", 
                                 nGroups, nSamps, nMol, nCores, B = NULL, 
                                 quantroAlpha = NULL, verbose = FALSE, 
                                 mua = NULL, siga = NULL, mub = NULL, sigb = NULL, 
                                 muOpt = NULL, sigOpt = NULL, muBG = NULL, 
                                 sigBG = NULL, muERR = NULL, sigERR = NULL)
{
    registerDoParallel(cores = nCores)
    workers <- getDoParWorkers()
    backend <- getDoParName()
    version <- getDoParVersion()
    
    N <- length(pDiffObject)
    lenP <- length(quantroAlpha)
    
    simResults <- foreach(i = 1:N) %dopar% {
        methTruth <- simulateMethTruth(nProbes = nMethProbesObject, 
                            nGroups = nGroups, pDiff = pDiffObject[i], 
                            pUp = pUp, verbose = verbose)
        if(nG==2){trueDiff = methTruth$methRange[,1] - methTruth$methRange[,2]}
        
        sim <- simulateMeth(methTruth, meth.platform = "methArrays", 
                            nSamps = nSamps, nMol = nMol, verbose = verbose, 
                            mua = mua, siga = siga, mub = mub, sigb = sigb, 
                            muOpt = muOpt, sigOpt = sigOpt, muBG = muBG, 
                            sigBG = sigBG, muERR = muERR, sigERR = sigERR)
         
        simDiffs = array(0, dim = c(nMethProbes, (2 + lenP)))

        Mset1 <- preprocessMeth(simulateMethObject = sim, method = "none")
        dmpPicksMset1 <- rowttests(as.matrix(Mset1$beta), sim$pd$Group)
        simDiffs[,1] <- dmpPicksMset1$dm
        
        Mset2 <- preprocessMeth(simulateMethObject = sim, method = "quantile")
        dmpPicksMset2 <- rowttests(as.matrix(Mset2$beta), sim$pd$Group)
        simDiffs[,2] <- dmpPicksMset2$dm
        
        simDiffs[,3:nMethods] <- sapply(1:lenP, function(x){
            Mset3 = preprocessMeth(simulateMethObject = sim, method = "quantro", 
                                   B = 100, alpha = quantroAlpha[x]) 
            dmpPicksMset3 = rowttests(as.matrix(Mset3$beta), sim$pd$Group)
            dmMset3 <- dmpPicksMset3$dm
        })
        
        return(cbind(trueDiff, simDiffs))
    }
    return(simResults)
}



calculateMSERandom <- function(simDiffsObject, nGroups = 2)
{
    nMethods <- dim(simDiffsObject)[2] - 1
    N <- dim(simDiffsObject)[3] 
    
    trueDiffObject <- sapply(1:N, function(x){ simDiffsObject[,1,x] })
    trueDiffObject <- abs(trueDiffObject)
    simDiffsObject <- abs(simDiffsObject)
    
    simMean = sapply(1:nMethods, function(x){ rowMeans(simDiffsObject[,(x+1),]) })
    simBias = sapply(1:nMethods, function(x){ rowMeans(simDiffsObject[,(x+1),] - trueDiffObject) })
    simMSE = sapply(1:nMethods, function(x){ rowMeans( (simDiffsObject[,(x+1),] - trueDiffObject )^2 ) })
    simVar <- colMeans(abs(simMSE)) - colMeans((simBias)^2)
    
    list(colMeans(abs(simMean)), colMeans(abs(simBias)), colMeans((simBias)^2), simVar, colMeans(simMSE))
}


calculateLossRandom <- function(simDiffsObject, nGroups = 2, topK = NULL)
{
    nMethods <- dim(simDiffsObject)[2] - 1
    N <- dim(simDiffsObject)[3]
    
    trueDiffObject <- sapply(1:N, function(x){ simDiffsObject[,1,x] })
    trueDiffObject <- abs(trueDiffObject)
    simDiffsObject <- abs(simDiffsObject)
    
    if(is.null(topK)){ stop("Must provide topK value.") }
    
    simLossMat = overlapMat <- array(0, dim = c(N, nMethods))
    for(i in 1:N){
        
        trueRank <- rank(-trueDiffObject[,i], ties.method = "first")
        
        methodsRank <- sapply(1:nMethods, function(x){
            rank(-simDiffsObject[,(x+1),i], ties.method = "first")
        })
        
        trueRankMat <- trueRank[order(trueRank)]
        methodsRankMat <- methodsRank[order(trueRank),]
        
        simLossMat[i,] <- sapply(1:nMethods, function(x){
            sum((trueRankMat[1:topK] - methodsRankMat[1:topK, x])^2)
        })    
        
        overlapMat[i,] <- sapply(1:nMethods, function(x){	
            length(intersect(trueRankMat[1:topK], methodsRankMat[1:topK, x]))
        })
    }	
    list("simLossMat" = simLossMat, "overlapMat" = overlapMat)
}


