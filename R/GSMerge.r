################## transform the result predlist to a matrix ###################
#' @export resultsMerge
resultsMerge <- function(predList){
  total_i <- dim.data.frame(predList)[2]
  predMatrix <- rbind(predList[[1]],predList[[2]])
  
  if (total_i > 2){ 
    for(i in 3:total_i){
      predMatrix <- rbind(predMatrix,predList[[i]])
    }
  }else{
    predMatrix <- predMatrix
  }
  predMatrix
}

############################# merge function ###########################
#' @export GSMerge
GSMerge <- function(predResMat, ratio, autoOptimize = F ){
  methodCount <- ncol(predResMat) - 1
  ## normalization
  rowCount <- nrow(predResMat)
  
  ## Auto optimize
  if(autoOptimize == T){
    
    # probMat <- transValMat2ProbMat(evalMat = predResMat,BestIndividuals = BestIndividuals)
    corMat <- cor(predResMat[,-1])
    evalresult <- G2PEvaluation(realScores = predResMat[,1],predScores = predResMat[,2:ncol(predResMat)],evalMethod = "pearson",topAlpha = 1:90)
    a <- evalresult$corMethods[1,]
    setMax <- which(a == max(a))
    b <- rank(a) - rank(corMat[,setMax])
    setPair <- which(b == max(b))[1]
    
    c <- apply(predResMat[,c((1+setMax),(1+setPair))],1,mean)

    mergePredRes <- as.matrix(c)
    mergePredRes <- cbind(predResMat[,c(1,setMax + 1,setPair + 1)],mergePredRes)
    colnames(mergePredRes)[ncol(mergePredRes)] <- "merge"
  }else{
    ## User optimize 
    mergePredRes <- as.matrix(predResMat[,-1] %*% diag(ratio))
    mergePredRes <- as.matrix(apply(mergePredRes,1,sum))/sum(ratio)
    colnames(mergePredRes) <- "merge"
    mergePredRes <- cbind(predResMat,mergePredRes)
  }
  mergePredRes
}


