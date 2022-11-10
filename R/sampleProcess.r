########### generated positive and negative samples for training ####################
#' @export sampleClassify
sampleClassify <- function(phenotype, posPercentage = 0.4, BestIndividuals = c("top", "middle", "buttom")){
  trainSampleSize <- length(phenotype)
  posSampleSize <- round( posPercentage*trainSampleSize )
  
  if( BestIndividuals == "top" ){
    trainPhenOrder <- order( phenotype, decreasing = TRUE )
  }else if( BestIndividuals == "buttom" ){
    trainPhenOrder <- order( phenotype, decreasing = FALSE )   
  }else if( BestIndividuals == "middle"){
    trainPhenOrder <- order( abs(phenotype), decreasing = FALSE ) 
  }
  posIdx <- trainPhenOrder[1:posSampleSize]
  negIdx <- setdiff( c(1:trainSampleSize), posIdx )
  
  res <- list( posSampleIndex = posIdx, negSampleIndex = negIdx )
  res
}
######################## generate train idx and test idx ##########################
#' @export cvSampleIndex
cvSampleIndex <- function( sampleNum, cross = 5, seed = 1,randomSeed = FALSE ) {
  if(randomSeed == TRUE){
    seed <- randomSeed()
  }
  cv <- cross
  resList <- list()
  
  # leave-one-out
  if( cv == sampleNum ){
    vec <- 1:sampleNum
    for( i in 1:sampleNum ){
      resList[[i]] <- list( trainIdx = vec[-i], testIdx = i, cvIdx = i)
    }
  }else {
    #random samples 
    set.seed(seed)
    index <- sample(1:sampleNum, sampleNum, replace = FALSE )
    step = floor( sampleNum/cv )
    
    start <- NULL
    end <- NULL
    train_sampleNums <- rep(0, cv)
    for( i in c(1:cv) ) {
      start <- step*(i-1) + 1
      end <- start + step - 1
      if( i == cv ) 
        end <- sampleNum
      
      testIdx <- index[start:end]
      trainIdx <- index[-c(start:end)]
      resList[[i]] <- list( trainIdx = trainIdx, testIdx = testIdx, cvIdx = i)
    }
  }
  names(resList) <- paste0("cv",1:cross)
  resList
}