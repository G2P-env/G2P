####################################################### Fit regression model ###################################################
#' @export fit.BGLR
#' @import BGLR
fit.BGLR <- function( trainMarkerMat, trainPheVal, predictMarkerMat, modelMethods,trainX = NULL, predX = NULL,  
                      outputModel = FALSE, nIter = 1500, burnIn = 500, thin = 5,verbose = FALSE,
                      saveAt = "", S0 = NULL, df0 =5, R2 = 0.5, weights = NULL,
                      rmExistingFiles = TRUE, groups=NULL){
  suppressMessages(require("BGLR"))
  if (!is.matrix(predictMarkerMat)) {
    predictMarkerMat <- t(as.matrix(predictMarkerMat))
  }
  numT <- nrow(predictMarkerMat)
  MarkerMat <- rbind(predictMarkerMat,trainMarkerMat)
  pheVal <- c(rep(NA,numT),trainPheVal)
  fix <- rbind(predX,trainX)
  
  checkModel <- modelMethods %in% c("BayesA", "BayesB", "BayesC", "BL", "BRR")
  if( ! checkModel ) {
    stop("Error: not defined model!")
  }
  
  if (!is.null(trainX) & !is.null(predX)) {
    ETA=list(list(X = MarkerMat, model= modelMethods),
             list(X = fix, model = "FIXED"))
  }else{
    ETA <- list(list(X=MarkerMat, model= modelMethods))
  }
  
  BGLRModel.fit <- BGLR( y = pheVal, ETA=ETA,  nIter = nIter, burnIn = burnIn, thin = thin, 
                         saveAt = saveAt, S0 = S0, df0 = df0, R2 = R2, weights = weights,
                         verbose = verbose, rmExistingFiles = rmExistingFiles, groups = groups )
  
  Res <- BGLRModel.fit$yHat[1:numT]
  if(outputModel){
    Res <- list(model = BGLRModel.fit,predictRes = Res)
  }
  Res
}
##
trainModel_RRBLUP <- function( markerMat, phenVec,X = NULL){
  phen_answer<-mixed.solve(phenVec, Z=markerMat, K=NULL, SE = FALSE, return.Hinv=FALSE,X = X)
  beta <- phen_answer$beta
  phD <- phen_answer$u
  e <- as.matrix(phD)
  return( list(beta = beta, e = e, phD = phD) )
}

### Fit Regression Model  
#' @import glmnet spls rrBLUP hglm hglm.data
#' @export GSReModel

GSReModel <- function(markers, pheVal, modelMethods,NAImpute = F,
                      K = 8, eta = 0.7, select = "pls2", fit = "simpls", scale.x = FALSE, scale.y = FALSE, eps = 1e-4, maxstep = 100,
                      alpha = 1,
                      X = NULL,family = gaussian(link = identity), lambda = NULL, tol.err = 1e-6, tol.conv = 1e-8, weights = NULL,
                      ...){  
  
  if (sum(is.na(pheVal)) != 0) {
    NAIdx <- which(is.na(pheVal) == T)
    if (NAImpute) {
      pheVal[NAIdx] <- mean(pheVal,na.rm =T)
    }else{
      pheVal <- pheVal[-NAIdx]
      markers <- markers[-NAIdx,]
      if (!is.null(X)) {
        X <- X[-NAIdx,]
      }
    }  
  }
  checkModel <- modelMethods %in% c("RRBLUP","LASSO","SPLS","bigRR")
  if( ! checkModel ) {
    stop("Error: not defined models for implementing GSReModel!")
  }
  #   if (modelMethods %in% c("BayesA", "BayesB", "BayesC", "BL", "BRR")){
  #     BGLRmethods <- modelMethods
  #     modelMethods <- "BGLRModel"
  #   }
  switch(modelMethods,
         #        BGLRModel  = trainedPredictModel_BGLR(trainMarkerMat = markers, trainedPhenVec = pheVal, modelMethods = BGLRmethods ,nIter = nIter, burnIn = burnIn, thin = thin, 
         #                                                saveAt = saveAt, S0 = S0, df0 = df0, R2 = R2, weights = weights,verbose = verbose, rmExistingFiles = rmExistingFiles, groups=groups),
         RRBLUP  = trainModel_RRBLUP(markerMat = markers, phenVec = pheVal,X = X),
         LASSO = trainModel_LASSO(markers,pheVal,alpha = alpha, ...),
         SPLS = trainModel_spls(markers,pheVal,K = K,eta = eta,select = select,fit = fit,scale.x = scale.x,scale.y = scale.y,eps = eps,maxstep = maxstep, ...),
         bigRR = trainModel_bigRR(markers = markers, X = X, pheVal = pheVal, weight = weights,family = family, lambda = lambda, tol.err = tol.err,tol.conv = tol.conv, ...)
  )
  
}

################################  LASSO ##############################
trainModel_LASSO <- function(markers,pheVal,alpha = 1, ...){
  suppressMessages(require("glmnet"))
  #glmnet fits a lasso model when we specify that alpha=1
  LASSO.fit <- glmnet(y=pheVal,x=markers,alpha=1, ...)
  #cv.glmnet finds the lambda value such that the the cvm value is the minimum
  cv <- cv.glmnet(y = pheVal, x=markers)
  LASSO_Res <- list(LASSO.fit = LASSO.fit,cv = cv)
  LASSO_Res
}

############################### spls #################################
trainModel_spls <- function(markers,pheVal,K = 8,eta = 0.7,select = "pls2",fit = "simpls",scale.x = FALSE,scale.y = FALSE,eps = 1e-4, maxstep = 100, ...){
  suppressMessages(require("spls"))
  f <- spls(markers,pheVal,K = K,eta = eta,select = select,fit = fit,scale.x =scale.x,scale.y = scale.y,eps = eps,trace = FALSE, ...)
  f
}

############################## bigRR #################################
trainModel_bigRR <- function (markers, X = NULL, pheVal, weight = NULL,
                              family = gaussian(link = identity), lambda = NULL, tol.err = 1e-6, 
                              tol.conv = 1e-8, ...){
  suppressMessages(require("hglm.data"))
  suppressMessages(require("hglm"))
  if(is.null(X)){
    X <- as.matrix(rep(1,length(pheVal)))
  }
  
  y <- pheVal
  Z <- markers
  
  Call <- match.call()
  if (!(is.matrix(X))) 
    stop("X should be a matrix.")
  if (!(is.matrix(Z))) 
    stop("Z should be a matrix.")
  if (!(is.vector(y))) 
    stop("y should be a vector.")
  
  if (any(is.na(y))) {
    naidx <- which(is.na(y))
    y <- y[-naidx]
    X <- X[-naidx,]
    Z <- Z[-naidx,]
  }
  
  N <- n <- nrow(X)
  p <- ncol(X)
  k <- ncol(Z)
  if (N != nrow(Z) | N != length(y)) 
    stop("Sizes of y, X, and Z are not all equal.")
  if (is.null(weight)) w <- rep(1, k) else w <- weight
  
  #G <- crossprod(sqrt(w)*t(Z)) ## bug fixed 111201 -- Xia
  wZt <- sqrt(w)*t(Z)
  G <- crossprod(wZt)
  ############ Bending to allow for p<n problems -- Lars (Xia added SVD)
  if (k < n) {
    eigen.values <- eigen(G)$values
    min.eigen <- min(eigen.values)
    if (min.eigen < tol.err) G <- G + diag(N)*(abs(min.eigen) + tol.err) 
  }
  ##
  #invG <- solve(G)
  #L <- t(chol(G))
  svdG <- svd(G)
  L <- svdG$u %*% diag(sqrt(svdG$d))
  invG <- tcrossprod(svdG$v %*% diag(1/svdG$d), svdG$u)
  phi0 <- sa0 <- 1
  if (is.null(lambda)) {
    hm <- hglm(y = y, X = X, Z = L, family = family, conv = tol.conv, bigRR = TRUE,...) ## checked with old emme code, conv = 1e-6 removed -- Xia
  }
  else {
    start.beta = c(rep(0, p))
    start.v = c(rep(0.01, n))
    start.lambda = lambda
    start.sigma2e = 1
    cat("Only 1 iteration applied for fixed lambda")
    hm <- hglm(y = y, X = X, Z = L, family = family, startval = c(start.beta, 
                                                                  start.v, start.lambda, start.sigma2e), maxit = 1,bigRR = TRUE, ...)
  }
  phi <- as.numeric(hm$varFix)
  sa <- as.numeric(hm$varRanef)
  a <- L%*%hm$ranef
  tZinvG <- crossprod(Z, invG)
  u <- (w*tZinvG)%*%a
  qu <- GCV <- NULL
  result <- list(phi = phi, lambda = hm$varRanef, beta = hm$fixef, hglm = hm,
                 u = u, leverage = qu, GCV = GCV, Call = Call, y = y, X = X)
  class(result) <- "bigRR"
  return(result)
}

############################### Fit a machine learning  model ############################
#' @export GSmachine
#' @import e1071 randomForest  
GSmachine <- function(markers, pheVal, modelMethods ="SVC", posPercentage = 0.4, BestIndividuals = c("top"),
                      ntree = 500,NAImpute = T,
                      nodesize = NULL, kernel = c("linear"), gamma = 1, cost = 2^(-9), ...){
  
  if (sum(is.na(pheVal)) != 0) {
    NAIdx <- which(is.na(pheVal) == T)
    if (NAImpute) {
      pheVal[NAIdx] <- mean(pheVal,na.rm =T)
    }else{
      pheVal <- pheVal[-NAIdx]
      markers <- markers[-NAIdx,]
    }  
  }
  
  if(is.null(nodesize)){
    if(modelMethods == "RFC"){
      nodesize <- 1
    }else if(modelMethods == "RFR"){
      nodesize <- 5
    }
  }else{
    nodesize <- nodesize
  }
  suppressMessages(require("e1071"))
  suppressMessages(require("randomForest")) 
  if( !modelMethods%in% c("SVR","SVC","RFR","RFC") ) {
    stop("Error: not defined category")
  }
  if (modelMethods %in%  c("SVC","RFC")){
    posNegSampleList <- sampleClassify(phenotype = pheVal ,posPercentage = posPercentage ,BestIndividuals = BestIndividuals )
    markers <- markers[c(posNegSampleList$posSampleIndex,posNegSampleList$negSampleIndex),]
    pheVal <-  as.factor( c( rep("1", length(posNegSampleList$posSampleIndex)), rep("0", length(posNegSampleList$negSampleIndex)) ) )
    
  }
  
  if(modelMethods %in%  c("SVR","SVC")){
    modelMethods <- "svm"
  }
  
  if(modelMethods %in% c("RFR","RFC")){
    modelMethods <- "RF"
  }
  switch(modelMethods,
         svm = svm(x= markers, y = pheVal,kernel = kernel,cost=cost,gamma = gamma,probability = TRUE, ...),
         RF = randomForest( x = markers, y = pheVal, ntree = ntree, importance = F,nodesize = nodesize), ...)
}

##############################################################################################################################3
########################### the prediction of genomic selection ###########################
#' @export predictGS
predictGS <- function(testMat, trainModel,predX = NULL,modelMethods = "SVC"){
  ########## check the methods of GS
  checkModel <- modelMethods %in% c("RRBLUP","SVR","SVC","RFR","RFC","LASSO","SPLS","bigRR")
  if( ! checkModel ) {
    stop("Error: not defined models for implementing GS Model")
  }
  
  Methods <- modelMethods
  
  ####### check testset 
  if(!is.matrix(testMat)) {
    testMat <- rbind(testMat,testMat)
    testnum <- 1 
  }else{
    testnum <- nrow(testMat)
  }
  
  #############
  if (!is.null(predX)) {
    if (modelMethods %in% c("RRBLUP","SVR","RFR","LASSO","SPLS","bigRR")){
      predresult <- switch(Methods,
                           #BGLRModel = {mkr.effs <- as.numeric(trainModel$ETA[[1]]$b); testMat %*% mkr.effs},
                           RRBLUP = {pred_phenVec <-testMat %*% trainModel$e;beta <- matrix(trainModel$beta,nrow = ncol(predX));beta <- predX %*% beta;as.numeric(pred_phenVec[,1]) + as.numeric(beta)},
                           SVR = {testMat <- cbind(testMat,predX);predict( trainModel, testMat)},
                           RFR = {testMat <- cbind(testMat,predX);predict( trainModel,testMat )},
                           LASSO = {testMat <- cbind(testMat,predX);pred_LASSO(trainModel,testMat)},
                           SPLS = {testMat <- cbind(testMat,predX);pred_SPLS(trainModel,testMat,type="fit")},
                           bigRR = {beta <- matrix(trainModel$beta,nrow = ncol(predX));beta <- predX %*% beta;as.numeric(trainModel$beta) + as.matrix(testMat %*% trainModel$u)}
      )
    }
    else if (modelMethods %in% c("SVC","RFC")){
      predresult <- switch(modelMethods,
                           SVC = {testMat <- cbind(testMat,predX);obj_pred <- predict(trainModel,testMat, probability = TRUE); as.matrix(attr(obj_pred, "probabilities")[,"1"])},
                           RFC = {testMat <- cbind(testMat,predX);predict(trainModel, testMat, type= "vote" )[,"1"]})
    }
  }else{
    if (modelMethods %in% c("RRBLUP","SVR","RFR","LASSO","SPLS","bigRR")){
      predresult <- switch(Methods,
                           #BGLRModel = {mkr.effs <- as.numeric(trainModel$ETA[[1]]$b); testMat %*% mkr.effs},
                           RRBLUP = {pred_phenVec <-  testMat %*% trainModel$e; as.numeric(pred_phenVec[,1]) + as.numeric(trainModel$beta)},
                           SVR = predict( trainModel, testMat),
                           RFR = predict( object = trainModel,testMat ),
                           LASSO = pred_LASSO(trainModel,testMat),
                           SPLS = pred_SPLS(trainModel,testMat,type= "fit"),
                           bigRR = as.numeric(trainModel$beta + testMat %*% trainModel$u)
      )
    }
    else if (modelMethods %in% c("SVC","RFC")){
      predresult <- switch(modelMethods,
                           SVC = {obj_pred <- predict(trainModel,testMat, probability = TRUE); as.matrix(attr(obj_pred, "probabilities")[,"1"])},
                           RFC = predict(trainModel, testMat, type= "vote" )[,"1"])
    }
  }
  predresult[1:testnum]
}


################### Modeling and predicting using RR ########################
#' @export fit.RR 
#' @import rrBLUP
fit.RR <- function(trainMarkerMat, trainPheVal, predictMarkerMat, NAImpute = F,trainX = NULL,predX = NULL,cpus = 1,outputModel = FALSE){
  suppressMessages(require("rrBLUP"))
  if (sum(is.na(trainPheVal)) != 0) {
    NAIdx <- which(is.na(trainPheVal) == T)
    if (NAImpute) {
      trainPheVal[NAIdx] <- mean(trainPheVal,na.rm =T)
    }else{
      trainPheVal <- trainPheVal[-NAIdx]
      trainMarkerMat <- trainMarkerMat[-NAIdx,]
      trainX <- trainX[-NAIdx,]
    }  
  }
  if (!is.matrix(predictMarkerMat)) {
    predictMarkerMat <- t(as.matrix(predictMarkerMat))
  }
  
  if(!is.null(trainX) & !is.null(predX)){
    yield_answer <- kinship.BLUP( y = trainPheVal, G.train = trainMarkerMat, G.pred= predictMarkerMat, K.method="RR",X = trainX, n.core = cpus)
    beta <- matrix(yield_answer$beta,nrow = ncol(trainX))
    beta <- predX %*% beta
    Res <-  as.numeric(yield_answer$g.pred) + as.numeric(beta)
  }else{
    yield_answer <- kinship.BLUP( y = trainPheVal, G.train = trainMarkerMat, G.pred= predictMarkerMat, K.method="RR", n.core = cpus)
    beta <- yield_answer$beta
    Res <-  as.numeric(yield_answer$g.pred) + as.numeric(beta)
  }
  
  if(outputModel){
    Res <- list(model = yield_answer,predictRes = Res)
  }
  Res
}

################### Modeling and predicting using BRNN#############################
#' @export fit.BRNN
#' @import brnn
fit.BRNN <- function(trainMarkerMat, trainPheVal, predictMarkerMat, NAImpute = F,outputModel = FALSE,verbose=TRUE, neurons=4, epochs=30, cpus = 1, ...){
  suppressMessages(require("brnn"))
  
  if (sum(is.na(trainPheVal)) != 0) {
    NAIdx <- which(is.na(trainPheVal) == T)
    if (NAImpute) {
      trainPheVal[NAIdx] <- mean(trainPheVal,na.rm =T)
    }else{
      trainPheVal <- trainPheVal[-NAIdx]
      trainMarkerMat <- trainMarkerMat[-NAIdx,]
    }  
  }
  
  if (!is.matrix(predictMarkerMat)) {
    predictMarkerMat <- t(as.matrix(predictMarkerMat))
  }
  X <- rbind(trainMarkerMat,predictMarkerMat)
  trainNum <- nrow(trainMarkerMat)
  yTRN <- trainPheVal
  n<-nrow(X) 
  p<-ncol(X)
  
  for(i in 1:ncol(X)){ (X[,i]<-X[,i]-mean(X[,i]))/sd(X[,i])}
  G<-tcrossprod(X)/ncol(X)
  trainMarkerMat <- G[1:trainNum,]
  predictMarkerMat <- G[-(1:trainNum),]
  
  NN<-brnn(y=yTRN,x=trainMarkerMat,neurons=neurons, epochs=epochs,verbose=verbose,cores = cpus,...)  
  Res <- predict(NN, newdata = predictMarkerMat)
  if(outputModel){
    Res <- list(model = NN,predictRes = Res)
  }
  Res
}
####################### LASSO pred ######################### 
pred_LASSO <- function(trainModel,testMat){
  predict(object = trainModel$LASSO.fit,testMat,s = trainModel$cv$lambda.min)
}

###################### pls and spls###########################
pred_SPLS <- function(trainModel,testMat,type = "fit"){
  predict( trainModel,testMat,type=type )
}

## RKHS #########
#' @export fit.RKHS
#' @import BGLR
fit.RKHS <- function(trainMarkerMat, trainPheVal, predictMarkerMat, NAImpute = F,trainX = NULL,predX = NULL,outputModel = FALSE,nIter = 1500, burnIn = 500, thin = 5,
                     saveAt = "", S0 = NULL, df0 =5, R2 = 0.5, weights = NULL,
                     verbose = FALSE, rmExistingFiles = TRUE, groups=NULL){
  suppressMessages(require("BGLR"))
  
  if (sum(is.na(trainPheVal)) != 0) {
    NAIdx <- which(is.na(trainPheVal) == T)
    if (NAImpute) {
      trainPheVal[NAIdx] <- mean(trainPheVal,na.rm =T)
    }else{
      trainPheVal <- trainPheVal[-NAIdx]
      trainMarkerMat <- trainMarkerMat[-NAIdx,]
      trainX <- trainX[-NAIdx,]
    }  
  }
  
  if (!is.matrix(predictMarkerMat)) {
    predictMarkerMat <- t(as.matrix(predictMarkerMat))
  }
  numT <- nrow(predictMarkerMat)
  MarkerMat <- rbind(predictMarkerMat,trainMarkerMat)
  fix <- rbind(predX,trainX)
  pheVal <- c(rep(NA,numT),trainPheVal)
  X <-scale(x = MarkerMat,center=TRUE,scale=TRUE)
  G <-tcrossprod(X)/ncol(X)
  
  if (!is.null(fix)) {
    ETA=list(list(K=G, model= "RKHS"),
             list(X=fix, model = "FIXED"))
  }else{
    ETA=list(list(K=G, model= "RKHS"))
  }
  RKHS.Fit <- BGLR(y=pheVal, ETA=ETA,  nIter = nIter, burnIn = burnIn, thin = thin, 
                   saveAt = saveAt, S0 = S0, df0 = df0, R2 = R2, weights = weights,
                   verbose = verbose, rmExistingFiles = rmExistingFiles, groups = groups)
  Res <- RKHS.Fit$yHat[1:numT]
  if(outputModel){
    Res <- list(model = RKHS.Fit,predictRes = Res)
  }
  Res
}