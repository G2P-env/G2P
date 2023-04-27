\name{GSmachine}
\alias{GSmachine}

\title{
Fit several machine learning models
}
\description{
The function fits several machine learning models of genomic selection such as SVR, SVC, RFR and RFC.
}
\usage{
GSmachine(markers, pheVal, modelMethods ="SVC", posPercentage = 0.4, 
          BestIndividuals = c("top"),ntree = 500,NAImpute = T,
          nodesize = NULL, kernel = c("linear"), gamma = 1, cost = 2^(-9), ...)
}

\arguments{
  \item{markers}{(numeric, matrix) row represents sample and column represents SNP (feature).Genotypes generally be coded as {0,1,2};0 represents AA(homozygote),2 represents BB(homozygote) and 1 represents AB(heterozygote); missing (NA) alleles are not allowed.}
  \item{pheVal}{(numeric) the phenotypic value of each sample.}
  \item{posPercentage}{(numeric,[0,1]) phenotype of extreme individuals which expected, default 0.4.}
  \item{BestIndividuals}{(character) the position of expected phenotype in whole phenotypic data set."top","buttom" or "middle",default "top".}
  \item{NAImpute}{(logical) if true, the missing values in phenotype will be imputed by mean otherwise be removed.}
  \item{modelMethods}{(character) alternative machine learning models. "SVR" and "SVC" from SVM, "RFR" and "RFC" from RF.}
  \item{ntree}{(integer) randomforest parameter. Number of trees to grow. This should not be set to too small a number, to ensure that every input row gets predicted at least a few times,default 500.}
  \item{nodesize}{(integer) randomforest parameter Minimum size of terminal nodes. Setting this number larger causes smaller trees to be grown (and thus take less time). Note that the default values are different for classification (1) and regression (5).}
  \item{kernel}{(numeric) svm parameter the kernel used in training and predicting. You might consider changing some of the following parameters, depending on the kernel type.(linear,polynomial,sigmoid,radial). Default "linear".}
  \item{gamma}{(numeric)svm parameter parameter needed for all kernels except linear, default 1.}
  \item{cost}{(numeric)svm parameter cost of constraints violation, default: 2^(-9), it is the 'C'-constant of the regularization term in the Lagrange formulation.}
}
\details{
SVM (support vector machine) and RF (random forest) models can be fitted in this function, including 2 classification (RFC, SVC) and 2 regression (SVR, SVC) models.
}
\value{
A fitted model
}
\references{
Meyer, D., Dimitriadou, E., Hornik, K., Weingessel, A., Leisch, F., Chang, C., et al. (2014). e1071: Misc functions of the Department of Statistics (e1071), TU Wien. R package version 1(3)

Liaw, A., and Wiener, M. (2002). Classification and regression by randomForest. R news 2(3), 18-22.

Classification and regression based on a forest of trees using random inputs, based on Breiman (2001) <DOI:10.1023/A:1010933404324>.
}

\seealso{
\pkg{e1017}\cr
\pkg{randomForest}
}

\examples{
## load example data ##
data("cubic")

## Fit RFR model ##
RFR_model <- GSmachine(markers = Markers[1:200,], pheVal = pheData$DTT[1:200], modelMethods = "RFR")

## Fit classification model (RFC) ##
RFC_model <- GSmachine(markers = Markers[1:200,], pheVal = pheData$DTT[1:200], modelMethods = "RFC",
                           posPercentage = 0.4, ntree = 500)
}