\name{fit.BRNN}
\alias{fit.BRNN}

\title{
Fitting BRNN model and performing prediction.
}

\description{
The function fits Bayesian regularization neural networks (BRNN) model and export the prediction value of testing sets. 
}
\usage{
fit.BRNN(trainMarkerMat, trainPheVal, predictMarkerMat, 
         NAImpute = FALSE,outputModel = FALSE,verbose=TRUE, 
         neurons=4, epochs=30,cpus = 1, ...)
}

\arguments{
  \item{trainMarkerMat}{(numeric, matrix) row represents sample and column represents SNP (feature).Genotypes generally be coded as {0,1,2};0 represents AA(homozygote),2 represents BB(homozygote) and 1 represents AB(heterozygote); missing (NA) alleles are not allowed.}
  \item{trainPheVal}{(numeric) the phenotypic value of each sample.}
  \item{predictMarkerMat}{(numeric, matrix) a new set of genotypic data that need to predict. row represents sample and column represents SNP (feature).Genotypes generally be coded as {0,1,2};0 represents AA(homozygote),2 represents BB(homozygote) and 1 represents AB(heterozygote), the coding pattern must be consistent with training data. missing (NA) alleles are not allowed.}
  \item{NAImpute}{(logical) if true, the missing values in phenotype will be imputed by mean otherwise be removed.}
  \item{cpus}{(integer) number of cpu cores to use for calculations (only available in UNIX-like operating systems), default 1.}
  \item{neurons}{(integer) indicates the number of neurons,defult 4.}
  \item{epochs}{(integer) maximum number of epochs(iterations) to train, default 30.}
  \item{verbose}{(logical) if TRUE, will print iteration history.}
  \item{outputModel}{(logical) if true, return the list of training model.}
}

\value{
a list including fitted model and prediction results (outputModel = TRUE).\cr
\item{model }{details of fitted model}
\item{predictRes }{prediction results of testing set}

a array of prediction results (outputModel = FALSE).
}
\references{
Pérez-Rodríguez, P., Gianola, D., Weigel, K., Rosa, G., and Crossa, J. (2013). Technical Note: An R package for fitting Bayesian regularized neural networks with applications in animal breeding. Journal of animal science 91(8), 3522-3531.
}

\seealso{
\pkg{brnn}
}
\examples{
## Load example data ##
data("cubic")

## use BRNN model to modeling and predict ##
RKHS_Res <- fit.RKHS(trainMarkerMat = Markers[1:200,], trainPheVal = pheData$DTT[1:200], 
                     predictMarkerMat = Markers[201:400,],outputModel = T,
                     nIter = 1500, burnIn = 500)
}