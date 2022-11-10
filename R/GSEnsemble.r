GSstacking <- function(predMat,evalMethods,by = 0.1,evaluation = T,topAlpha = 15, ...){
  globalMethods <-  c("pearson", "kendall", "spearman", "MSE","R2")
  if (evalMethods %in% globalMethods) {
    A <- "corMethods"
  }else{
    A <- evalMethods
  }
  rowNum <- nrow(predMat)
  colNum <- ncol(predMat)
  observeVal <- predMat[,1]
  
  subEnsemble <- predMat[,c(1,2,3)]
  
  radioEnsem <- GSMerge(predResMat = subEnsemble,ratio = c(0,1),autoOptimize = F)
  
  cycle <- seq(0,1,by = by)
  cycleLength <- length(cycle)
  
  for(ii in cycle[-1]){
    ensembleRes <- GSMerge(predResMat = subEnsemble,ratio = c(ii,1-ii),autoOptimize = F)
    radioEnsem <- cbind(radioEnsem,as.matrix(ensembleRes[,4]))
  }
  
  
  evalGYSSRatio <- G2PEvaluation(realScores = radioEnsem[,1], predScores = radioEnsem[,4:(cycleLength + 3)], 
                              evalMethod = evalMethods,topAlpha = topAlpha, ...)
  if(evalMethods == "MSE"){
    maxEnsembleIdex <- order(evalGYSSRatio[[A]])[1]  
  }else{
    maxEnsembleIdex <- order(evalGYSSRatio[[A]],decreasing = T)[1]
  }
  
  maxEnsemble <- evalGYSSRatio[[A]][maxEnsembleIdex]
  Index <- which(evalGYSSRatio[[A]] == maxEnsemble)
  if (length(Index) > 1 ) {
    Index <- Index[1]
  }
  weightmodel1 <- (Index - 1)/(cycleLength-1)
  weightmodel2 <- (cycleLength - Index)/(cycleLength-1)
  ensembleRes <- radioEnsem[,Index + 3]
  weight <- c(weightmodel1,weightmodel2)
  
  for(i in 4:colNum){
    subEnsemble <- predMat[,1:i]
    radioEnsem <- GSMerge(predResMat = subEnsemble,ratio = c(0*weight,1),autoOptimize = F)
    for(ii in cycle[-1]){
      ensembleRes <- GSMerge(predResMat = subEnsemble,ratio = c(ii*weight,1-ii),autoOptimize = F)
      radioEnsem <- cbind(radioEnsem,as.matrix(ensembleRes[,i+1]))
    }
    evalGYSSRatio <- G2PEvaluation(realScores = radioEnsem[,1], predScores = radioEnsem[,(i+1):(i+cycleLength)], 
                                evalMethod = evalMethods,topAlpha = topAlpha, ...)
    if(evalMethods == "MSE"){
      maxEnsembleIdex <- order(evalGYSSRatio[[A]])[1]  
    }else{
      maxEnsembleIdex <- order(evalGYSSRatio[[A]],decreasing = T)[1]
    }
    
    maxEnsemble <- evalGYSSRatio[[A]][maxEnsembleIdex]
    Index <- which(evalGYSSRatio[[A]] == maxEnsemble)[1]
    if (length(Index) > 1 ) {
      Index <- Index[1]
    }
    weightmodel1 <- (Index - 1)/(cycleLength-1)
    weightmodel2 <- (cycleLength - Index)/(cycleLength-1)
    ensembleRes <- radioEnsem[,Index + i]
    
    weight <- c(weight*weightmodel1,weightmodel2)
  }
  names(weight) <- colnames(predMat)[-1]
  finalMat <- cbind(predMat,as.matrix(ensembleRes))
  colnames(finalMat)[colNum+1] <- "Ensemble"
  
  if(evaluation){
    evalRes <- G2PEvaluation(realScores = finalMat[,1], predScores = finalMat[,2:(colNum+1)], 
                          evalMethod = evalMethods,topAlpha = topAlpha, ... )
    res <- list(Weights = weight,finalMat = finalMat,evalRes = evalRes)
  }else{
    res <- list(Weights = weight,finalMat = finalMat)
  }
  res
}

############################# Ensemble function ###########################
#' @export GSEnsemble
GSEnsemble <- function(predMat, nrandom = 10, evalMethods, by = 0.1, evaluation = T, topAlpha = 15, ...){
  rowNum <- nrow(predMat)
  colNum <- ncol(predMat)
  methodsNames <- colnames(predMat)[-1]
  methodsNames1 <- c(colnames(predMat)[-1],"Ensemble")
  
  weightMat <- matrix(NA,colNum-1,nrandom)
  evalMat <-  matrix(NA,colNum,nrandom)
  rownames(weightMat) <- methodsNames
  rownames(evalMat) <- methodsNames1
  
  for(i in 1:nrandom){
    a <- sample(1:(colNum-1),(colNum-1),replace = F)
    newRes <- predMat[,2:colNum]
    newRes <- newRes[,a]
    newRes <- cbind(as.matrix(predMat[,1]),newRes)
    colnames(newRes)[1] <- colnames(predMat)[1]
    
    res <- GSstacking(predMat = newRes, evalMethods = evalMethods,by = by,evaluation = evaluation,topAlpha = topAlpha)
    if (evalMethods %in% c("pearson", "kendall", "spearman", "MSE","R2")) {
      Methods <- "corMethods"  
    }else{
      Methods <- evalMethods
    }
    weightMat[,i] <- res$Weights[methodsNames]
    evalMat[,i] <- res$evalRes[[Methods]][,methodsNames1]
  }
#   weightMat <- weightMat[,-1]
#   pearsonMat <- pearsonMat[,-1]
  colnames(weightMat) <- paste0("rep",1:nrandom)
  colnames(evalMat) <- paste0("rep",1:nrandom)
  
  index <- which(evalMat["Ensemble",] == max(evalMat["Ensemble",]))
  if(length(index) > 1){
    index <- index[1]
  }
  bestWeight <- weightMat[,index]
  finalMat <- GSMerge(predResMat = predMat,ratio = weightMat[,index],autoOptimize = F)
  colnames(finalMat)[colNum + 1] <- "Ensemble"
  
  if(evaluation){
    evalRes <- G2PEvaluation(realScores = finalMat[,1], predScores = finalMat[,2:(colNum+1)], 
                          evalMethod = evalMethods,topAlpha = topAlpha, ...)
    res <- list(BestWeight = bestWeight,finalMat = finalMat,evalRes = evalRes,weightMat = weightMat,evalMat = evalMat)
  }else{
    res <- list(BestWeight = bestWeight,finalMat = finalMat,weightMat = weightMat,evalMat = evalMat)
  }
  res
}

