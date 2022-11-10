G2Pnofix <- function(trainMarker,trainPheno,testMarker,testPheno = NULL,modelMethods ="BayesA",outputModel = FALSE,  # main parameters
                nIter = 1500, burnIn = 500, thin = 5, 
                saveAt = "", S0 = NULL, df0 =5, R2 = 0.5, weights = NULL,
                verbose = FALSE, rmExistingFiles = TRUE, groups=NULL,importance = FALSE,    # # # BGLR method parameters
                posPercentage = 0.4,BestIndividuals = c("top"),ntree = 500,nodesize = NULL,kernel = c("linear"),gamma = 1, cost = 2^(-9),  # machine learing parameters
                K = 8,eta = 0.7,select = "pls2",fit = "simpls",scale.x = FALSE,scale.y = FALSE,eps = 1e-4,maxstep = 100, # SPLS parameters
                alpha = 1,X = NULL,family = gaussian(link = identity), lambda = NULL, tol.err = 1e-6, tol.conv = 1e-8,
                epochs = 30, neurons = 4,
                ...){   # LASSO parameters
  if(is.null(testPheno)){
    ntestSample <- nrow(testMarker)
    if (is.null(ntestSample)) {
      testPheno <- rep(0,1)
    }else{
      testPheno <- rep(0,ntestSample)
    }
  }
  allModelMethod <- c('SPLS','LASSO','SVC','SVR','RFC','RFR','RRBLUP','bigRR','BayesA','BayesB','BayesC','BL','BRR','RKHS','RR','BRNN')
  modelMethods <- modelMethods
  ## judge 
  if (!all(modelMethods %in% allModelMethod)) {
    stop("The provided list of modelMethods is out of G2P model list, please check!")
  }
  
  require("pbapply")
  pboptions(type = "timer")
  result <- pblapply(modelMethods, function(ii){cat(ii,"is modeling ...","\n");
    switch(ii,
           RR = fit.RR(trainMarker,trainPheno,testMarker,outputModel = TRUE),
           BRNN = fit.BRNN(trainMarker,trainPheno,testMarker,outputModel = TRUE,verbose = verbose, neurons = neurons, epochs = epochs, ...),
           RKHS = fit.RKHS(trainMarker,trainPheno,testMarker,outputModel = TRUE, nIter = nIter, burnIn = burnIn,thin = thin,
                           saveAt = saveAt, S0 = S0, df0 = df0, R2 = R2, weights = weights,verbose = verbose, rmExistingFiles = rmExistingFiles, groups=groups),
           BayesA = fit.BGLR(trainMarker, trainPheno, testMarker, modelMethods = "BayesA",outputModel = TRUE,nIter = nIter, burnIn = burnIn, thin = thin,
                             saveAt = saveAt, S0 = S0, df0 =df0, R2 = R2, weights = weights,verbose = verbose, rmExistingFiles = rmExistingFiles, groups=groups),
           BayesB = fit.BGLR(trainMarker, trainPheno, testMarker, modelMethods = "BayesB",outputModel = TRUE,nIter = nIter, burnIn = burnIn, thin = thin,
                             saveAt = saveAt, S0 = S0, df0 =df0, R2 = R2, weights = weights,verbose = verbose, rmExistingFiles = rmExistingFiles, groups=groups),
           BayesC = fit.BGLR(trainMarker, trainPheno, testMarker, modelMethods = "BayesC",outputModel = TRUE,nIter = nIter, burnIn = burnIn, thin = thin,
                             saveAt = saveAt, S0 = S0, df0 =df0, R2 = R2, weights = weights,verbose = verbose, rmExistingFiles = rmExistingFiles, groups=groups),
           BL = fit.BGLR(trainMarker, trainPheno, testMarker, modelMethods = "BL",outputModel = TRUE,nIter = nIter, burnIn = burnIn, thin = thin,
                         saveAt = saveAt, S0 = S0, df0 =df0, R2 = R2, weights = weights,verbose = verbose, rmExistingFiles = rmExistingFiles, groups=groups),
           BRR = fit.BGLR(trainMarker, trainPheno, testMarker, modelMethods = "BRR",outputModel = TRUE,nIter = nIter, burnIn = burnIn, thin = thin,
                          saveAt = saveAt, S0 = S0, df0 =df0, R2 = R2, weights = weights,verbose = verbose, rmExistingFiles = rmExistingFiles, groups=groups),
           RRBLUP = {resModel <- GSReModel(modelMethods = "RRBLUP",markers = trainMarker,pheVal = trainPheno,
                                           K = K,eta = eta,select = select,fit = fit,scale.x = scale.x,scale.y = scale.y,eps = eps,
                                           alpha = alpha,X = X,family = family, lambda = lambda, tol.err = tol.err, tol.conv = tol.conv,maxstep = maxstep);
           predictRes <- predictGS(testMat = testMarker,trainModel = resModel,modelMethods = "RRBLUP");
           list(model = resModel,predictRes = predictRes)},
           LASSO = {resModel <- GSReModel(modelMethods = "LASSO",markers = trainMarker,pheVal = trainPheno,
                                          K = K,eta = eta,select = select,fit = fit,scale.x = scale.x,scale.y = scale.y,eps = eps,
                                          alpha = alpha,X = X,family = family, lambda = lambda, tol.err = tol.err, tol.conv = tol.conv,maxstep = maxstep);
           predictRes <- predictGS(testMat = testMarker,trainModel = resModel,modelMethods = "LASSO");
           list(model = resModel,predictRes = predictRes)},
           SPLS = {resModel <- GSReModel(modelMethods = "SPLS",markers = trainMarker,pheVal = trainPheno,
                                         K = K,eta = eta,select = select,fit = fit,scale.x = scale.x,scale.y = scale.y,eps = eps,
                                         alpha = alpha,X = X,family = family, lambda = lambda, tol.err = tol.err, tol.conv = tol.conv,maxstep = maxstep);
           predictRes <- predictGS(testMat = testMarker,trainModel = resModel,modelMethods = "SPLS");
           list(model = resModel,predictRes = predictRes)},
           bigRR = {resModel <- GSReModel(modelMethods = "bigRR",markers = trainMarker,pheVal = trainPheno,
                                          K = K,eta = eta,select = select,fit = fit,scale.x = scale.x,scale.y = scale.y,eps = eps,
                                          alpha = alpha,X = X,family = family, lambda = lambda, tol.err = tol.err, tol.conv = tol.conv,maxstep = maxstep);
           predictRes <- predictGS(testMat = testMarker,trainModel = resModel,modelMethods = "bigRR");
           list(model = resModel,predictRes = predictRes)},
           SVC = {resModel <- GSmachine(markers = trainMarker,pheVal = trainPheno,posPercentage = posPercentage ,BestIndividuals = BestIndividuals,
                                        modelMethods = "SVC" ,ntree = ntree ,nodesize = nodesize,
                                        kernel = kernel,gamma = gamma, cost = cost);
           predictRes <- predictGS(testMat = testMarker,trainModel = resModel,modelMethods = "SVC");
           list(model = resModel,predictRes = predictRes)},
           
           RFC = {resModel <- GSmachine(markers = trainMarker,pheVal = trainPheno,posPercentage = posPercentage ,BestIndividuals = BestIndividuals,
                                        modelMethods = "RFC" ,ntree = ntree ,nodesize = nodesize,
                                        kernel = kernel,gamma = gamma, cost = cost);
           predictRes <- predictGS(testMat = testMarker,trainModel = resModel,modelMethods = "RFC");
           list(model = resModel,predictRes = predictRes)},
           
           SVR = {resModel <- GSmachine(markers = trainMarker,pheVal = trainPheno,posPercentage = posPercentage ,BestIndividuals = BestIndividuals,
                                        modelMethods = "SVR" ,ntree = ntree ,nodesize = nodesize,
                                        kernel = kernel,gamma = gamma, cost = cost);
           predictRes <- predictGS(testMat = testMarker,trainModel = resModel,modelMethods = "SVR");
           list(model = resModel,predictRes = predictRes)},
           
           RFR = {resModel <- GSmachine(markers = trainMarker,pheVal = trainPheno,posPercentage = posPercentage ,BestIndividuals = BestIndividuals,
                                        modelMethods = "RFR" ,ntree = ntree ,nodesize = nodesize,
                                        kernel = kernel,gamma = gamma, cost = cost);
           predictRes <- predictGS(testMat = testMarker,trainModel = resModel,modelMethods = "RFR");
           list(model = resModel,predictRes = predictRes)}
    )
  })
  
  ModelList <- lapply(result, function(x){a <- x[[1]]})
  PredResList <- lapply(result, function(x){b <- x[[2]]})
  lengtList <- length(result)
  realPheno <- as.matrix(testPheno)
  resMat <- as.matrix(realPheno)
  for(i in 1:lengtList){
    resMat <- cbind(resMat,as.matrix(PredResList[[i]]))
  }
  colnames(resMat) <- c("realPhenScore",modelMethods)
  rownames(resMat) <- rownames(testMarker)
  names(ModelList) <- modelMethods
  finalRes <- list(model = ModelList, predictResMat = resMat)
  # predscores <- predscores[,c("realPhenScore",modelMethods)]
  ####### output result
  if(outputModel){
    return(finalRes)
  }
  else{
    return(resMat)
  }
}


## 
## complex Genotype to Phenotype (CG2P) with fix ########################################
#' @export G2P
G2P <- function(markers,data,fix = NULL,trait,trainIdx,predIdx,modelMethods ="BayesA",outputModel = FALSE,NAImpute = FALSE,  # main parameters
                nIter = 1500, burnIn = 500, thin = 5, 
                saveAt = "", S0 = NULL, df0 =5, R2 = 0.5, weights = NULL,
                verbose = FALSE, rmExistingFiles = TRUE, groups=NULL,importance = FALSE,    # # # BGLR method parameters
                posPercentage = 0.4,BestIndividuals = c("top"),ntree = 500,nodesize = NULL,kernel = c("linear"),gamma = 1, cost = 2^(-9),  # machine learing parameters
                K = 8,eta = 0.7,select = "pls2",fit = "simpls",scale.x = FALSE,scale.y = FALSE,eps = 1e-4,maxstep = 100, # SPLS parameters
                alpha = 1,family = gaussian(link = identity), lambda = NULL, tol.err = 1e-6, tol.conv = 1e-8,
                epochs = 30, neurons = 4,
                ...){   # LASSO parameters
  
  trainMarker <- markers[trainIdx,]
  testMarker <- markers[predIdx,]
  trainPheno <- data[trainIdx,trait]
  testPheno <- data[predIdx,trait]
  
  if (NAImpute) {
    trainPheno[which(is.na(trainPheno) == T)] <- mean(trainPheno,na.rm = T) 
  }
  
  if(!is.null(fix)){
    ## treatment of fix predictor
    fix.data.mat <- eval(parse(text = paste0("model.matrix(~",paste0(fix,collapse = " + "),"-1,data = data)")))
    
    trainX <- fix.data.mat[trainIdx,]
    predX <- fix.data.mat[predIdx,]
    
    allModelMethod <-  c('SPLS','LASSO','SVC','SVR','RFC','RFR','RRBLUP','bigRR','BayesA','BayesB','BayesC','BL','BRR','RKHS','RR','BRNN')
    
    ## judge 
    if (!all(modelMethods %in% allModelMethod)) {
      stop("The provided list of modelMethods is out of G2P model list, please check!")
    }
    
    require("pbapply")
    pboptions(type = "timer")
    result <- pblapply(modelMethods, function(ii){cat(ii,"is modeling ...","\n");
      switch(ii,
             RR = fit.RR(trainMarker,trainPheno,testMarker,trainX = trainX,predX = predX,outputModel = TRUE),
             BRNN = {trainMarker <- cbind(trainMarker,trainX);testMarker <- cbind(testMarker,predX);fit.BRNN(trainMarker,trainPheno,testMarker,outputModel = TRUE,verbose = verbose, neurons = neurons, epochs = epochs, ...)},
             RKHS = fit.RKHS(trainMarker,trainPheno,testMarker,outputModel = TRUE,trainX = trainX,predX = predX,nIter = nIter, burnIn = burnIn,thin = thin,
                             saveAt = saveAt, S0 = S0, df0 = df0, R2 = R2, weights = weights,verbose = verbose, rmExistingFiles = rmExistingFiles, groups=groups),
             BayesA = fit.BGLR(trainMarker, trainPheno, testMarker, modelMethods = "BayesA",trainX,predX,outputModel = TRUE,nIter = nIter, burnIn = burnIn, thin = thin,
                               saveAt = saveAt, S0 = S0, df0 =df0, R2 = R2, weights = weights,verbose = verbose, rmExistingFiles = rmExistingFiles, groups=groups),
             BayesB = fit.BGLR(trainMarker, trainPheno, testMarker, modelMethods = "BayesB",trainX,predX,outputModel = TRUE,nIter = nIter, burnIn = burnIn, thin = thin,
                               saveAt = saveAt, S0 = S0, df0 =df0, R2 = R2, weights = weights,verbose = verbose, rmExistingFiles = rmExistingFiles, groups=groups),
             BayesC = fit.BGLR(trainMarker, trainPheno, testMarker, modelMethods = "BayesC",trainX,predX,outputModel = TRUE,nIter = nIter, burnIn = burnIn, thin = thin,
                               saveAt = saveAt, S0 = S0, df0 =df0, R2 = R2, weights = weights,verbose = verbose, rmExistingFiles = rmExistingFiles, groups=groups),
             BL = fit.BGLR(trainMarker, trainPheno, testMarker, modelMethods = "BL",trainX,predX,outputModel = TRUE,nIter = nIter, burnIn = burnIn, thin = thin,
                           saveAt = saveAt, S0 = S0, df0 =df0, R2 = R2, weights = weights,verbose = verbose, rmExistingFiles = rmExistingFiles, groups=groups),
             BRR = fit.BGLR(trainMarker, trainPheno, testMarker, modelMethods = "BRR",trainX,predX,outputModel = TRUE,nIter = nIter, burnIn = burnIn, thin = thin,
                            saveAt = saveAt, S0 = S0, df0 =df0, R2 = R2, weights = weights,verbose = verbose, rmExistingFiles = rmExistingFiles, groups=groups),
             RRBLUP = {resModel <- GSReModel(modelMethods = "RRBLUP",markers = trainMarker,pheVal = trainPheno,
                                             K = K,eta = eta,select = select,fit = fit,scale.x = scale.x,scale.y = scale.y,eps = eps,
                                             alpha = alpha,X = trainX,family = family, lambda = lambda, tol.err = tol.err, tol.conv = tol.conv,maxstep = maxstep);
             predictRes <- predictGS(testMat = testMarker,trainModel = resModel,modelMethods = "RRBLUP",predX = predX);
             list(model = resModel,predictRes = predictRes)},
             LASSO = {markers <- cbind(markers,fix.data.mat);trainMarker <- markers[trainIdx,];resModel <- GSReModel(modelMethods = "LASSO",markers = trainMarker,pheVal = trainPheno,
                                                                                                                     K = K,eta = eta,select = select,fit = fit,scale.x = scale.x,scale.y = scale.y,eps = eps,
                                                                                                                     alpha = alpha,family = family, lambda = lambda, tol.err = tol.err, tol.conv = tol.conv,maxstep = maxstep);
             predictRes <- predictGS(testMat = testMarker,trainModel = resModel,modelMethods = "LASSO",predX = predX);
             list(model = resModel,predictRes = predictRes)},
             SPLS = {markers <- cbind(markers,fix.data.mat);trainMarker <- markers[trainIdx,];resModel <- GSReModel(modelMethods = "SPLS",markers = trainMarker,pheVal = trainPheno,
                                                                                                                    K = K,eta = eta,select = select,fit = fit,scale.x = scale.x,scale.y = scale.y,eps = eps,
                                                                                                                    alpha = alpha,family = family, lambda = lambda, tol.err = tol.err, tol.conv = tol.conv,maxstep = maxstep);
             predictRes <- predictGS(testMat = testMarker,trainModel = resModel,modelMethods = "SPLS",predX = predX);
             list(model = resModel,predictRes = predictRes)},
             bigRR = {resModel <- GSReModel(modelMethods = "bigRR",markers = trainMarker,pheVal = trainPheno,
                                            K = K,eta = eta,select = select,fit = fit,scale.x = scale.x,scale.y = scale.y,eps = eps,
                                            alpha = alpha,X = trainX,family = family, lambda = lambda, tol.err = tol.err, tol.conv = tol.conv,maxstep = maxstep);
             predictRes <- predictGS(testMat = testMarker,trainModel = resModel,modelMethods = "bigRR",predX = predX);
             list(model = resModel,predictRes = predictRes)},
             SVC = {markers <- cbind(markers,fix.data.mat);trainMarker <- markers[trainIdx,];resModel <- GSmachine(markers = trainMarker,pheVal = trainPheno,posPercentage = posPercentage ,BestIndividuals = BestIndividuals,
                                                                                                                   modelMethods = "SVC" ,ntree = ntree ,nodesize = nodesize,
                                                                                                                   kernel = kernel,gamma = gamma, cost = cost);
             predictRes <- predictGS(testMat = testMarker,trainModel = resModel,modelMethods = "SVC",predX = predX);
             list(model = resModel,predictRes = predictRes)},
             
             RFC = {markers <- cbind(markers,fix.data.mat);trainMarker <- markers[trainIdx,];resModel <- GSmachine(markers = trainMarker,pheVal = trainPheno,posPercentage = posPercentage ,BestIndividuals = BestIndividuals,
                                                                                                                   modelMethods = "RFC" ,ntree = ntree ,nodesize = nodesize,
                                                                                                                   kernel = kernel,gamma = gamma, cost = cost);
             predictRes <- predictGS(testMat = testMarker,trainModel = resModel,modelMethods = "RFC",predX = predX);
             list(model = resModel,predictRes = predictRes)},
             
             SVR = {markers <- cbind(markers,fix.data.mat);trainMarker <- markers[trainIdx,];resModel <- GSmachine(markers = trainMarker,pheVal = trainPheno,posPercentage = posPercentage ,BestIndividuals = BestIndividuals,
                                                                                                                   modelMethods = "SVR" ,ntree = ntree ,nodesize = nodesize,
                                                                                                                   kernel = kernel,gamma = gamma, cost = cost);
             predictRes <- predictGS(testMat = testMarker,trainModel = resModel,modelMethods = "SVR",predX = predX);
             list(model = resModel,predictRes = predictRes)},
             
             RFR = {markers <- cbind(markers,fix.data.mat);trainMarker <- markers[trainIdx,];resModel <- GSmachine(markers = trainMarker,pheVal = trainPheno,posPercentage = posPercentage ,BestIndividuals = BestIndividuals,
                                                                                                                   modelMethods = "RFR" ,ntree = ntree ,nodesize = nodesize,
                                                                                                                   kernel = kernel,gamma = gamma, cost = cost);
             predictRes <- predictGS(testMat = testMarker,trainModel = resModel,modelMethods = "RFR",predX = predX);
             list(model = resModel,predictRes = predictRes)}
      )
    })
    
    ModelList <- lapply(result, function(x){a <- x[[1]]})
    PredResList <- lapply(result, function(x){b <- x[[2]]})
    lengtList <- length(result)
    realPheno <- as.matrix(testPheno)
    resMat <- as.matrix(realPheno)
    for(i in 1:lengtList){
      resMat <- cbind(resMat,as.matrix(PredResList[[i]]))
    }
    colnames(resMat) <- c("realPhenScore",modelMethods)
    rownames(resMat) <- rownames(testMarker)
    names(ModelList) <- modelMethods
    # predscores <- predscores[,c("realPhenScore",modelMethods)]
    ####### output result
    if(outputModel){
      finalRes <- list(model = ModelList, predictResMat = resMat)
    }else{
      finalRes <- resMat
    }
    finalRes
  }else{
    finalRes <- G2Pnofix(trainMarker,trainPheno,testMarker,testPheno = testPheno,modelMethods = modelMethods,outputModel = outputModel,
                         nIter = nIter, burnIn = burnIn, thin = thin, NAImpute = NAImpute,
                         saveAt = saveAt, S0 = S0, df0 =df0, R2 = R2, weights = weights,
                         verbose = verbose, rmExistingFiles = rmExistingFiles, groups=groups,importance = importance, 
                         posPercentage = posPercentage,BestIndividuals =BestIndividuals,ntree = ntree,nodesize = nodesize,kernel = kernel,gamma = gamma, cost =cost, 
                         K = K,eta = eta,select = select,fit = fit,scale.x = scale.x,scale.y = scale.y,eps = eps,maxstep = maxstep, 
                         alpha = alpha,family = family, lambda = lambda, tol.err = tol.err, tol.conv = tol.conv,
                         epochs = epochs, neurons = neurons,
                         ...)
  }
  return(finalRes)
}