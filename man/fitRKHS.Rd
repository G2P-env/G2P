\name{fitRKHS}
\alias{fitRKHS}

\title{
Fitting RKHS model and performing prediction.
}

\description{
The function fits Reproducing kernel Hilbert space (RKHS) model and export the prediction value of testing sets. 
}
\usage{
fit.RKHS(trainMarkerMat, trainPheVal, predictMarkerMat, 
         NAImpute = F,trainX = NULL,predX = NULL,outputModel = FALSE,
         nIter = 1500, burnIn = 500, thin = 5,
         saveAt = "", S0 = NULL, df0 =5, R2 = 0.5, weights = NULL,
         verbose = FALSE, rmExistingFiles = TRUE, groups=NULL)
}

\arguments{
  \item{trainMarkerMat}{(numeric, matrix) genotypic data of training set, row represents sample and column represents SNP (feature).Genotypes generally be coded as {0,1,2};0 represents AA(homozygote),2 represents BB(homozygote) and 1 represents AB(heterozygote); missing (NA) alleles are not allowed.}
  \item{trainPheVal}{(numeric) the phenotypic value of each sample.}
  \item{outputModel}{(logical) if TRUE, return the list of trained model and prediction results, default FALSE.}
  \item{predictMarkerMat}{(numeric, matrix) genotypic data of testing set, row represents sample and column represents SNP (feature).Genotypes generally be coded as {0,1,2};0 represents AA(homozygote),2 represents BB(homozygote) and 1 represents AB(heterozygote); missing (NA) alleles are not allowed.}
  \item{NAImpute}{(logical) if true, the missing values in phenotype will be imputed by mean otherwise be removed.}
  \item{trainX}{additional data inputs of training sets as fixed effects instead of SNPs inputs only.}
  \item{predX}{additional data inputs of testing sets as fixed effects instead of SNPs inputs only.}
  \item{nIter,burnIn,thin}{(integer) the number of iterations, burn-in and thinning,default nIter 7000,burnIn 500,thin 5.}
  \item{saveAt}{(string) this may include a path and a pre-fix that will be added to the name of the files that are saved as the program runs,default "".}  
  \item{S0,df0}{(numeric) the scale parameter for the scaled inverse-chi squared prior assigned to the residual variance, only used with Gaussian outcomes.Default S0 NULL,df0 = 5.}
  \item{R2}{(numeric, (0,1)) the proportion of variance that one expects, a priori, to be explained by the regression. Only used if the hyper-parameters are not specified.Defult 0.5}
  \item{weights}{(numeric) a vector of weights, may be NULL. If weights is not NULL, the residual variance of each data-point is set to be proportional to the square of the weight. Only used with Gaussian outcomes.}
  \item{verbose}{(logical) if TRUE the iteration history is printed, default FALSE.}
  \item{rmExistingFiles}{(logical) if TRUE, removes existing output files from previous runs, default TRUE.}
  \item{groups}{(factor) a vector of the same length of y that associates observations with groups, each group will have an associated variance component for the error term.}
}

\value{
a list including fitted model and prediction results (outputModel = TRUE).\cr
\item{model }{details of fitted model}
\item{predictRes }{prediction results of testing set}

a array of prediction results (outputModel = FALSE).
}
\references{
Pérez, P., & de los Campos, G. (2014). BGLR: a statistical package for whole genome regression and prediction. Genetics, 198(2), 483-495.
}

\seealso{
\pkg{BGLR}
}

\examples{
## Load example data ##
data("cubic")

## use RKHS model to modeling and predict ##
RKHS_Res <- fit.RKHS(trainMarkerMat = Markers[1:200,], trainPheVal = pheData$DTT[1:200], 
                     predictMarkerMat = Markers[201:400,],nIter = 1500, burnIn = 500)
}