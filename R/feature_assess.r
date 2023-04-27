############################## Feature Selection ############################
#'@import rrBLUP randomForest
#'@export feature_assess 
feature_assess <- function(markers, phenotype, method = c("RRBLUP", "Gini", "Accuracy"), ntree = 500, importance = TRUE, posPercentage = 0.40, BestIndividuals = c("top", "middle", "buttom")){
  require("rrBLUP")
  require("randomForest")
  
  if(length(method) > 1){
    method = "rrBLUP"
  }
  if(length(BestIndividuals > 1)){
    BestIndividuals <- "top"
  }
  posNegSampleList <- sampleClassify(phenotype = phenotype, posPercentage = posPercentage, BestIndividuals = BestIndividuals )
  Xtrain <- markers[c(posNegSampleList$posSampleIndex,posNegSampleList$negSampleIndex),]
  YClasif <-  as.factor( c( rep("1", length(posNegSampleList$posSampleIndex)), rep("0", length(posNegSampleList$negSampleIndex)) ) )
  switch(method,
         RRBLUP = trainModel_RRBLUP(markers , phenotype)$phD,
         Accuracy =  randomForest(x = Xtrain, y = YClasif, ntree = ntree, importance = importance)$importance[,"MeanDecreaseAccuracy"],
         Gini = randomForest( x = Xtrain, y = YClasif, ntree = ntree, importance = importance)$importance[,"MeanDecreaseGini"]
  )
}
