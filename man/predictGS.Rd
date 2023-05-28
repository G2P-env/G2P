\name{predictGS}
\alias{predictGS}

\title{
Perform prediciton with fitted Models.
}
\description{
The function perform prediction from new inputs with fitted models, including RRBLUP, SVR, SVC, RFR, RFC, LASSO, SPLS and bigRR.
}
\usage{
predictGS(testMat, trainModel,predX = NULL,modelMethods = "SVC")
}

\arguments{
  \item{testMat}{(numeric, matrix) a new set of genotypic data that need to predict. row represents sample and column represents SNP (feature).Genotypes generally be coded as {0,1,2};0 represents AA(homozygote),2 represents BB(homozygote) and 1 represents AB(heterozygote), the coding pattern must be consistent with training data. missing (NA) alleles are not allowed.}
  \item{trainModel}{(model) the fitted model. The object must be consistent with modelMethods.}
  \item{predX}{(numeric, matrix) a design matrix as fixed effects, the formation must be consistent with fitted model used.}
  \item{modelMethods}{(character) the name of fitted model, including RRBLUP, SVR, SVC, RFR, RFC, LASSO, SPLS and bigRR.}
}

\value{
a list:\cr  
The prediction results of test set which predicted through already models.
}

\examples{
## Load example data ##
data("cubic")

## Fit RRBLUP model ##
RRBLUP_model <- GSReModel(markers = Markers[1:200,], pheVal = pheData$DTT[1:200], modelMethods = "RRBLUP")

## Predict 201-400 subset of all example data with already RRBLUP model ## 
res <- predictGS(testMat = Markers[201:400,], trainModel = RRBLUP_model, modelMethods = "RRBLUP")
}