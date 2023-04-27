#####################################  cross validation for complex genomic selection with fix ####################################
#' @export G2P
#' @import parallel pbapply
G2PCrossValidation <-function(cvSampleList = NULL,cross = 10,times = 1,seed = 1,cpus = 1, markers, data, trait,fix = NULL,
                               modelMethods ="SVC", outputModel = FALSE,NAImpute = FALSE,
                               nIter = 1500, burnIn = 500, thin = 5, 
                               saveAt = "", S0 = NULL, df0 =5, R2 = 0.5, weights = NULL,
                               verbose = FALSE, rmExistingFiles = TRUE, groups=NULL, importance = FALSE,
                               posPercentage = 0.4, BestIndividuals = c("top"), ntree = 500, nodesize = NULL, kernel = c("linear"), gamma = 1, cost = 2^(-9),
                               K = 8, eta = 0.7, select = "pls2", fit = "simpls", scale.x = FALSE, scale.y = FALSE, eps = 1e-4, trace = FALSE, maxstep = 100,  # SPLS parameters
                               alpha = 1,family = gaussian(link = identity), lambda = NULL, tol.err = 1e-6, tol.conv = 1e-8,
                               epochs = 30, neurons = 4,
                               ...){
  suppressMessages(require("parallel"))
  suppressMessages(require("pbapply"))
  
  sampleNum <- nrow(markers)
  phenotypeNum <- dim(data)[1]
  cvSampleList = cvSampleList;
  pheVal <- data[,trait]
  
  data = data;trait = trait;fix = fix;NAImpute = NAImpute;
  cross = cross; seed = seed; cpus = cpus; markers = markers;pheVal = pheVal;
  modelMethods = modelMethods; BestIndividuals = BestIndividuals;
  nIter = nIter; burnIn = burnIn; thin = thin; saveAt = saveAt; S0 = S0; df0 =df0; R2 = R2; weights = weights;posPercentage = posPercentage;
  verbose = verbose; rmExistingFiles = rmExistingFiles; groups=groups;ntree = ntree  ;nodesize = nodesize ;importance = importance;
  kernel = kernel;gamma = gamma; cost = cost;outputModel = outputModel;
  K = K;eta = eta;select = select;fit = fit;scale.x = scale.x;scale.y = scale.y;eps = eps;trace = trace;
  alpha = alpha;family = family; lambda = lambda; tol.err = tol.err; tol.conv = tol.conv;maxstep = maxstep;
  epochs = epochs; neurons = neurons;
  
  if(sampleNum != phenotypeNum) {
    stop("Marker count is not equal to phenotype count!")
  }
  
  if(!is.numeric(markers) | !is.numeric(pheVal)){
    stop("Marker or phenotype is not numeric, please check it!")  
  }
  
  
  cl <- makeForkCluster(cpus)
  cat(cpus," cores were used for cross validation ... \n")
  cat("Start cross validation ... \n")
  
  if(!is.null(cvSampleList)){
    cvSampleList <- cvSampleList
  }else{
    if (times == 1) {
      cvSampleList <- cvSampleIndex(sampleNum = sampleNum,cross = cross,seed = seed,randomSeed = F)
    }else if (times > 1){
      cvSampleListToatal <- c()
      for (i in 1:times) {
        cvSampleList <- cvSampleIndex(sampleNum = sampleNum,cross = cross,seed = seed,randomSeed = T)
        names(cvSampleList) <- paste0("Times_",i,"_",names(cvSampleList))
        cvSampleListToatal <- c(cvSampleListToatal,cvSampleList)
        cvSampleListToatal
      }
      cvSampleList <- cvSampleListToatal
    }
  }
  
  pboptions(type = "timer")
  results <- pblapply(1:length(cvSampleList),function(x){
    # library("BGLR")
    suppressMessages(library("G2P"))
    suppressMessages(library("brnn"))
    suppressMessages(library("glmnet"))
    suppressMessages(library("spls"))
    suppressMessages(library("pls"))
    suppressMessages(library("e1071"))
    suppressMessages(library("BGLR"))
    suppressMessages(library("rrBLUP"))
    suppressMessages(library("randomForest"))
    suppressMessages(library("hglm"))
    suppressMessages(library("hglm.data"))
    # source("./CG2P.r")
    
    cat("All needed package have loaded in all cores! \n")
    
    markers <- markers
    data <- data
    trait <- trait
    fix <- fix
    trainIdx <- cvSampleList[[x]]$trainIdx
    predIdx <-  cvSampleList[[x]]$testIdx 
    
    G2P(data = data,markers = markers,fix = fix,trait = trait,NAImpute = NAImpute,trainIdx = trainIdx,predIdx = predIdx,modelMethods = modelMethods ,BestIndividuals = BestIndividuals,
         nIter = nIter, burnIn = burnIn, thin = thin,  saveAt = saveAt, S0 = S0, df0 =df0, R2 = R2, weights = weights,posPercentage = posPercentage,
         verbose = verbose, rmExistingFiles = rmExistingFiles, groups=groups,ntree = ntree  ,nodesize = nodesize ,importance = importance,
         kernel = kernel,gamma = gamma, cost = cost,outputModel = outputModel,
         K = K,eta = eta,select = select,fit = fit,scale.x = scale.x,scale.y = scale.y,eps = eps,trace = trace,
         alpha = alpha,family = family, lambda = lambda, tol.err = tol.err, tol.conv = tol.conv,maxstep = maxstep,
         epochs = epochs, neurons = neurons,
         ...)
  },cl=cl) # lapply
  stopCluster(cl)
  cat(times," times",cross," fold cross validation is done! \n")
  if (times == 1) {
    names(results) <- paste0("CV",1:cross)
  }else{
    final_res <- list()
    length(final_res) <- times
    names(final_res) <- paste0("Rep",1:times)
    for (i in 1:times) {
      final_res[[i]] <- results[((i-1)*cross):(i*cross)]
      names(final_res[[i]]) <- paste0("CV",1:cross)
    }
    results <- final_res
  }
  results
}
