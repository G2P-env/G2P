################################################ eval ###############################################
meanNDCG <- function( realScores, predScores, topK = 10 ){
  resVec <- rep(0, topK )
  for( idx in 1:topK ){
    resVec[idx] <- NDCG( realScores, predScores, topK = idx )
  }
  meanNDCG <- mean(resVec)
  names(meanNDCG ) <- paste0("meanNDCG","_top",topK)
  return (meanNDCG)
}
##from plos one, 2015, A Ranking Approach to Genomic Selection
NDCG <- function( realScores, predScores, topK = 10){
  
  if( length(realScores) != length(predScores)){
    stop("Error: different length between realScores and predScores")
  }
  if( length(realScores) < topK ) {
    stop("Error: too large topK")
  }
  
  scoreMat <- cbind(realScores,predScores)
  scoreMatSortbyPred <- scoreMat[order(scoreMat[,2],decreasing = TRUE),]
  scoreMatSortByReal <- scoreMat[order(scoreMat[,1],decreasing = TRUE),]
  
  DCG <- rep(0, topK)
  IDCG <- rep(0, topK)
  for(idx in 1:topK){
    DCG[idx] <-  scoreMatSortbyPred[idx,1]/log(idx+1,2)
    IDCG[idx] <- scoreMatSortByReal[idx,1]/log(idx+1,2)
  }
  
  NDCG <- sum(DCG)/sum(IDCG) 
  names(NDCG) <- paste0("NDCG","_top",topK)
  return(NDCG)
}
##evaluation method: pearson, spearman, kendall, MSE
corEvaluation <- function( realScores, predScores, method = c("pearson", "kendall", "spearman", "MSE","R2"),BestIndividuals = "top",probIndex = NULL){
  # Probability handle
  if (!is.null(probIndex)) {
    if(BestIndividuals == "top"){
      realScores <- realScores
      predScores[,probIndex] <- predScores[,probIndex]
    }else if(BestIndividuals == "buttom"){
      realScores <- realScores
      predScores[,probIndex] <- 1 - predScores[,probIndex]
    }else if(BestIndividuals == "middle"){
      realScores <- abs(realScores)
      predScores[,probIndex] <- 1 - predScores[,probIndex] 
    }
  }else{
    realScores <- realScores
    predScores <- predScores
  }
  
  if( length(method) > 1 ){
    method <- method[1]
  }
  checkMethodType <- method %in% c("pearson", "kendall", "spearman", "MSE","R2")
  if( !checkMethodType ){
    stop("Error: undefined method in corEvaluation")
  }
  
  realScores <- as.matrix(realScores)
  predScores <- as.matrix(predScores)
  res <- ""
  if( (method == "pearson") | (method == "kendall") | (method == "spearman") ){
    res <- cor( realScores, predScores,  use="complete", method = method  )
  }else if(method == "MSE") {
    res <- apply(predScores,2,function(ii){
      deltaVec <- abs( realScores - ii)
      deltaVec <- deltaVec^2
      mean(deltaVec)})
  }else if(method == "R2"){
    res <-apply(predScores,2,function(ii){
      R2 <- summary(lm(realScores ~ ii))$r.squared })
  }
  
  res <- matrix(res,nrow = 1,dimnames = list(method,colnames(predScores)))
  res
}
##
## evaluation method alpha : pearson, spearman, kendall, MSE
corEvaluationAlpha <- function( realScores, predScores, method = c("pearson", "kendall", "spearman", "MSE","R2"),topAlpha = 15,BestIndividuals = "top",probIndex = NULL){
  topNum <- 1: round(length(realScores)*(topAlpha/100))
  if( length(method) > 1 ){
    method <- method[1]
  }
  checkMethodType <- method %in% c("pearson", "kendall", "spearman", "MSE","R2")
  if( !checkMethodType ){
    stop("Error: undefined method in corEvaluation")
  }
  realScores <- as.matrix(realScores)
  predScores <- as.matrix(predScores)
  mat <- cbind(realScores,predScores)
  
  if(BestIndividuals == "top"){
    mat <- mat[order(realScores,decreasing = T),]
  }else if(BestIndividuals == "buttom"){
    mat <- mat[order(realScores,decreasing = F),]
  }else if(BestIndividuals == "middle"){
    mat <- mat[order(abs(realScores),decreasing = F),]
  }
  
  mat <- mat[topNum,]
  
  realScores <- as.matrix(mat[,1])
  predScores <- as.matrix(mat[,-1])
  
  res <- corEvaluation(realScores, predScores, method = method,BestIndividuals = BestIndividuals,probIndex = probIndex)
  res
}

multiAlphaCor <- function(realScores, predScores, method = c("pearson", "kendall", "spearman", "MSE","R2"),topAlpha = 1:15,BestIndividuals = "top",probIndex = NULL){
  res <- lapply(method,function(meth){
    # cormat <- apply(predScores,2,function(x){
    cormat <- sapply(topAlpha,function(ii){
      corEvaluationAlpha(realScores = realScores,predScores = predScores,method = meth,topAlpha = ii,BestIndividuals = BestIndividuals, probIndex = probIndex)})
    cormat <- matrix(cormat,ncol = length(topAlpha))
    cormat <- t(cormat)
    rownames(cormat) <- paste0("top",topAlpha)
    colnames(cormat) <- colnames(as.matrix(predScores))
    cormat
  })
  names(res) <- method
  res
}
####
## function: RE and kappa from Heredity paper
## type:regression
## rank: top:(TRUE,TRUE); middle and down :(FALSE,FALSE)
## type: classfication
## rank: top:(TRUE,TRUE): middle and down :top:(FALSE,TRUE)
classEvaluation <- function(realScores, predScores, topAlpha = 15, Beta = 1,Probability = FALSE,evalMethod = "AUC",BestIndividuals = c("top", "middle", "buttom") ) {
  
  #if( is.null( names(realScores)) ){
  #names(realScores) <- paste("sample", 1:length(realScores), sep = "")
  #}
  if( length(realScores) != length(predScores) ){
    stop("Error: different length between realScores and predScores")
  }
  
  if( length(BestIndividuals) > 1 ) {
    BestIndividuals <- BestIndividuals[1]
  }
  
  realScores <- as.numeric(realScores)
  predScores <- as.numeric(predScores)
  total <- length(realScores)
  topK <- round( total*topAlpha/100 )
  classVec <- c( rep(1, topK), rep(0, total-topK))
  ##
  if(BestIndividuals == "top"){
    decreaseReal <- TRUE
    decreasePred <- TRUE
  }else if(BestIndividuals == "buttom"){
    if(Probability){
      decreaseReal <- FALSE
      decreasePred <- TRUE
    }else{
      decreaseReal <- FALSE
      decreasePred <- FALSE
    }
  }else if(BestIndividuals == "middle"){
    realScores <- abs(realScores)
    predScores <- abs(predScores)
    if(Probability){
      decreaseReal <- FALSE
      decreasePred <- TRUE
    }else{
      decreaseReal <- FALSE
      decreasePred <- FALSE
    }
  }
  
  scoreMat <- cbind( realScores, predScores )
  newScoreMat <- scoreMat[order(scoreMat[,1], decreasing = decreaseReal),] 
  newScoreMat <- cbind( newScoreMat, classVec )
  topRealMean <- mean( newScoreMat[1:topK,1] )
  #
  threshold <- newScoreMat[topK,1]
  #
  
  ##### RE ,kappa
  newScoreMat <- newScoreMat[order(newScoreMat[,2], decreasing = decreasePred),]
  # classVec <- c(rep(1,length(which(newScoreMat[,2] <= threshold))),rep(0, total-(length(which(newScoreMat[,2] <= threshold)))))
  newScoreMat <- cbind( newScoreMat, classVec )
  colnames(newScoreMat) <- c("real", "predicted", "realClass", "predClass")
  PredRealMean <- mean( newScoreMat[1:topK,1] ) 
  
  TP <- sum(newScoreMat[1:topK, 3])
  FP <- topK - TP
  FN <- topK - TP
  TN <- total - topK - FN
  Po <- (TP+TN)/total
  Pe <- ((FP+TN)/total)*((FN+TN)/total) + ((TP+FN)/total)*((TP+FP)/total)
  allRealMean <- mean( scoreMat[,1] )
  precision = TP/(TP + FP)
  recall = TP/(TP +FN)
  ###  the area under the receiver operating characteristics curve(AUC),and  the area under the precision-recall  curve (AUCpr)
  result <- sapply(evalMethod,function(one_method){
    switch(one_method,
           Kappa =  (Po-Pe)/(1-Pe),
           RE = ( PredRealMean - allRealMean)/(topRealMean - allRealMean),
           #AUC = roc.curve(scores.class0 = newScoreMat[newScoreMat[,3] == 1,2],scores.class1 = newScoreMat[newScoreMat[,3] == 0,2],curve = TRUE)$AUC ,
           AUC = roc(newScoreMat[,3],newScoreMat[,4],quiet = T)$auc[[1]],
           AUCpr = pr.curve(scores.class0 = newScoreMat[newScoreMat[,3] == 1,2],scores.class1 = newScoreMat[newScoreMat[,3] == 0,2],curve = TRUE,sorted = TRUE)$auc.integral,
           accuracy = (TP + TN)/total,
           F1 = (1 + Beta^2)*precision *recall/(precision + recall )
    )
  })
  names(result) <- paste(evalMethod,"_top", topAlpha, sep = "" )
  return(result) 
}
##################################### to evaluate for multiple parameters set
multiParameters <- function(realScores,predScores, Probability = TRUE,evalMethod = c("RE"),topAlpha ,Beta = 1,BestIndividuals = c("top"),probIndex = NULL){
  
  ############ calculate one or multiple topAlpha 
  multiTopAlpha <- function(realScores,predScores, Probability ,evalMethod ,topAlpha,Beta,BestIndividuals){
    ######################### one or multiple topAlpha for classifiction evaluation
    if (length(intersect(evalMethod,c("RE", "Kappa", "AUC","AUCpr","accuracy" ,"precision","recall","F1" ))) != 0){
      result <-  sapply(topAlpha,function(ii)
        classEvaluation(realScores = realScores,predScores = predScores,Probability = Probability,
                        evalMethod = evalMethod,topAlpha = ii,Beta = Beta,BestIndividuals = BestIndividuals)
      ) }
    
    ######################### one or multiple topAlpha for NDCG evaluation
    ################ set format of output
    if(length(evalMethod) > 1){
      result <- t(result)
    }else{
      result <- as.matrix(result)
    }
    dimnames(result) <- list(paste0("top",topAlpha),evalMethod )
    return(result)
  }
  
  
  predScores <- as.matrix(predScores)
  evalNum <- 1:ncol(predScores)
  
  ## class ##
  if(!is.null(probIndex)){
    classMethSite <- probIndex
    evalNum1 <- evalNum[-classMethSite]
    multiPredresultS <- lapply(evalNum1,function(ii){
      multiTopAlpha(realScores = realScores,predScores = predScores[,ii],Probability = FALSE,evalMethod = evalMethod ,topAlpha = topAlpha,Beta = Beta,BestIndividuals = BestIndividuals )
    })
    classPredresult <- lapply(classMethSite,function(ii){
      multiTopAlpha(realScores = realScores,predScores = predScores[,ii],Probability = TRUE,evalMethod = evalMethod ,topAlpha = topAlpha,Beta = Beta,BestIndividuals = BestIndividuals )
    })
    multiPredresult <- list()
    length(multiPredresult) <- ncol(predScores)
    multiPredresult[evalNum1] <- multiPredresultS
    multiPredresult[classMethSite] <- classPredresult
  }else{
    ############ evaluate one or multiple prediction result
    multiPredresult <- lapply(evalNum,function(ii){
      multiTopAlpha(realScores = realScores,predScores = predScores[,ii],Probability = FALSE,evalMethod = evalMethod ,topAlpha = topAlpha,Beta = Beta,BestIndividuals = BestIndividuals )
    })
  }
  ######### set format of output
  result <- lapply(evalMethod,function(ii){
    one_evalMethod <- sapply(evalNum,function(jj) multiPredresult[[jj]][,ii])
    if(length(topAlpha) > 1){
      colnames(one_evalMethod) <- colnames(predScores)
    }else{
      one_evalMethod <- matrix(one_evalMethod,nrow = 1,dimnames = list(paste0("top",topAlpha),colnames(predScores)))
    }
    one_evalMethod
  })
  names(result) <- evalMethod
  return(result)
}
################################## evaluateGS ###################################
#' @export G2PEvaluation
#' @import pROC PRROC
G2PEvaluation <- function(realScores, predScores, evalMethod = "RE", Beta = 1, 
                          BestIndividuals = "top", topAlpha = 1:90, 
                          globalAlpha = F, probIndex = NULL){
  suppressMessages(require(pROC))
  suppressMessages(require(PRROC))
  ## data process
  ## remove NAs
  evalMat <- cbind(as.matrix(realScores),as.matrix(predScores))
  evalMat <- na.omit(evalMat)
  ## warning prediction results not converging
  check <- apply(as.matrix(evalMat[,2:ncol(evalMat)]),2,function(x){
    if (length(unique(x)) == 1) {
      1
    }else{
      0
    }
  })
  
  if (1 %in% check) {
    warning(paste0("The prediction resluts of model ",paste0(names(check)[check == 1]),collapse = ",")," is abnormal. Please check.")
    abn.idx <- which(check == 1)
    if (!is.null(probIndex)){
      for(i in 1:length(probIndex)){
        diff.idx <- probIndex[i] - abn.idx
        ## minus
        probIndex[i] <- probIndex[i] - length(diff.idx[diff.idx > 0])
      }
    }
    evalMat <- evalMat[,-(which(check == 1)+1)]
    if (!is.matrix(evalMat)) {
      stop("There is no predictions value and cannot be assessed!")
    }
  }
  
  # probMat <- transValMat2ProbMat(evalMat = evalMat,BestIndividuals = BestIndividuals)
  realScores <- evalMat[,1]
  predScores <- evalMat[,-1]
  ## Evaluation
  selectFun <- NULL
  predScores <- as.matrix(predScores)
  globalMethods <-  c("pearson", "kendall", "spearman", "MSE","R2")
  thresholdMethods <- c("RE", "Kappa", "AUC","AUCpr","accuracy" ,"precision","recall","F1" )
  NDCGMethods <- c("meanNDCG", "NDCG")
  ###################### correlation methods
  if (length(intersect(evalMethod,globalMethods)) != 0){
    if(globalAlpha){
      selectFun <- c(selectFun,"globalEvaluation","globalEvaluationAlpha")
      globalMethods  <- intersect(evalMethod,globalMethods)
    }else{
      selectFun <- c(selectFun,"globalEvaluation")
      globalMethods  <- intersect(evalMethod,globalMethods)
    }
  }
  
  ################# classification evaluation
  if (length(intersect(evalMethod,thresholdMethods)) != 0){
    selectFun <- c(selectFun,"thresholdEvaluation")
    thresholdMethods <- intersect(evalMethod,thresholdMethods)
  }
  
  result <- lapply(selectFun,function(one_fun){
    switch(one_fun,
           globalEvaluation = {corResult <- sapply(globalMethods,function(one_methods){corEvaluation( realScores, predScores, method = one_methods,BestIndividuals = BestIndividuals,probIndex = probIndex)});
           matrix(t(corResult),ncol= ncol(predScores),dimnames = list(globalMethods,colnames(predScores)))},
           globalEvaluationAlpha = multiAlphaCor(realScores,predScores,topAlpha = topAlpha,method = globalMethods,BestIndividuals = BestIndividuals,probIndex = probIndex),
           thresholdEvaluation = multiParameters(realScores, predScores, topAlpha = topAlpha, Probability = Probability,evalMethod = thresholdMethods,Beta =Beta , BestIndividuals = BestIndividuals, probIndex = probIndex)
    )})
  ############ the format of output
  finalresult <- list()
  nList <- length(result)
  nList <- c("one","two","three")[nList]
  if(length(result) != 0){
    #     id <- 1
    #     if(!is.list(result[[1]])){finalresult[["corMethods"]] <- result[[1]];id <- 2};finalresult <- c(finalresult,result[[id]])
    switch(nList,
           one = {if(!is.list(result[[1]])){finalresult[["corMethods"]] <- result[[1]]}else{finalresult <- c(finalresult,result[[1]])}},
           two = {if(!is.list(result[[1]])){finalresult[["corMethods"]] <- result[[1]];finalresult <- c(finalresult,result[[2]])}else{finalresult <- c(finalresult,result[[1]],result[[2]])}},
           three = {finalresult[["corMethods"]] <- result[[1]];finalresult <- c(finalresult,result[[2]],result[[3]])})
  }else{finalresult <- list()}
  #################### NDCG evaluation
  if (length(intersect(evalMethod,NDCGMethods)) != 0){
    selectFun <- intersect(evalMethod,NDCGMethods)
    result2 <- lapply(selectFun,function(one_fun){
      switch (one_fun,
              NDCG = apply(predScores,2,function(x){NDCGEvaluation(method = "NDCG",topAlpha = topAlpha,realScores = realScores,predScores = x)}),
              meanNDCG = apply(predScores,2,function(x){NDCGEvaluation(method = "meanNDCG",topAlpha = topAlpha,realScores = realScores,predScores = x)})
      )})
    names(result2) <- intersect(evalMethod,NDCGMethods)
    finalresult <- c(finalresult,result2)
  }
  ### RFC SVC
  # if ((!is.null(probIndex)) & ("corMethods" %in% names(finalresult))){
  #   probColIndex <- probIndex
  #   predScores <- as.matrix(predScores[,probColIndex])
  #   corResult <- sapply(globalMethods,function(one_methods){corEvaluation( realScores, predScores, method = one_methods,BestIndividuals = BestIndividuals,Probability = Probability)})
  #   corResult <- t(corResult)
  #   finalresult$corMethods[,probColIndex] <- corResult
  # }
  finalresult
}
NDCGEvaluation <- function(method = NULL,topAlpha,realScores,predScores, allNDCG = T){
  
  if(allNDCG == T){
    ndcg <- 1 : round(length(realScores)*(max(topAlpha)/100))
  }else{
    ndcg <- round(length(realScores)*(topAlpha/100))
  }
  
  if(method == "NDCG"){
    result <- sapply(ndcg,function(ii){
      NDCG(realScores = realScores,predScores = predScores,topK = ii)
      
    })
    if(allNDCG == F){
      names(result) <- paste0("NDCG_top",topAlpha)
    }
  }else if(method == "meanNDCG"){
    ndcg <- 1 : round(length(realScores)*(max(topAlpha)/100))
    NDCGEval <- sapply(ndcg,function(ii){
      NDCG(realScores = realScores,predScores = predScores,topK = ii)
    })
    result <- sapply(topAlpha,function(ii){
      iii <- round(ii*length(realScores)/100)
      sum(NDCGEval[1:iii])/iii
    })
    names(result) <- paste0("meanNDCG_top",topAlpha)
  }
  result
}