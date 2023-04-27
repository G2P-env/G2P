\name{GSReModel}
\alias{GSReModel}

\title{
Fit several common regression models
}

\description{
The function can fit several regression models of genomic selection including LASSO, RRBLUP, SPLS and bigRR.
}
\usage{
GSReModel(markers, pheVal, modelMethods,NAImpute = F,
          K = 8, eta = 0.7, select = "pls2", fit = "simpls",
          scale.x = FALSE, scale.y = FALSE, eps = 1e-4, maxstep = 100,
          alpha = 1,X = NULL,family = gaussian(link = identity), lambda = NULL,
          tol.err = 1e-6, tol.conv = 1e-8, weights = NULL, ...)
}

\arguments{
  \item{markers}{(numeric, matrix) row represents sample and column represents SNP (feature).Genotypes generally be coded as {0,1,2};0 represents AA(homozygote),2 represents BB(homozygote) and 1 represents AB(heterozygote); missing (NA) alleles are not allowed.}
  \item{pheVal}{(numeric) the phenotypic value of each sample.}
  \item{modelMethods}{(character) alternative regression models. "RRBLUP", "LASSO", "SPLS" and "bigRR".}
  \item{NAImpute}{(logical) if true, the missing values in phenotype will be imputed by mean otherwise be removed.}
  \item{K}{(integer) SPLS parameter: number of hidden components.}
  \item{eta}{(numeric) SPLS parameter: thresholding parameter. eta should be between 0 and 1.}
  \item{select}{(character) SPLS parameter: PLS algorithm for variable selection. Alternatives are "pls2" or "simpls". Default is "pls2".} 
  \item{fit}{(character) SPLS parameter: PLS algorithm for model fitting. Alternatives are "kernelpls", "widekernelpls", "simpls", or "oscorespls". Default is "simpls".}
  \item{scale.x}{(character) SPLS parameter: scale predictors by dividing each predictor variable by its sample standard deviation?}
  \item{scale.y}{(character) SPLS parameter: scale responses by dividing each response variable by its sample standard deviation?}
  \item{eps}{(character) SPLS parameter: an effective zero. Default is 1e-4.}
  \item{maxstep}{(integer) SPLS parameter: maximum number of iterations when fitting direction vectors. Default is 100.}
  \item{alpha}{(numeric) LASSO model parameter: the elasticnet mixing parameter.Detail in glmnet.}
  \item{X}{(numeric, matrix) design matrix related to the parameters not to be shrunk (i.e. fixed effects in the mixed model framework),defult no shrink.}
  \item{family}{the distribution family of y, see help('family') for more details. }
  \item{lambda}{the shrinkage parameter determines the amount of shrinkage. Default is NULL meaning that it is to be estimated along with other model parameters.}
  \item{tol.err}{internal tolerance level for extremely small values; default value is 1e-6.}
  \item{tol.conv}{tolerance level in convergence; default value is 1e-8.}
}
\details{
SPLS (Sparse partial least squares)\cr
RRBLUP (Ridge regression best linear unbiased prediction)\cr
LASSO (Least absolute shrinkage and selection operator)\cr
bigRR (big ridge regression)
}
\value{
A list model which is enable to predict.
}
\references{
Endelman, J.B. 2011. Ridge regression and other kernels for genomic selection with R package rrBLUP. Plant Genome 4:250-255. doi: 10.3835/plantgenome2011.08.0024

Chung, D., Chun, H., and Keles, S. (2012). Spls: sparse partial least squares (SPLS) regression and classification. R package, version 2, 1-1.

Friedman, J., Hastie, T., and Tibshirani, R. (2010). Regularization paths for generalized linear models via coordinate descent. Journal of statistical software 33(1), 1

hglm: A package for fitting hierarchical generalized linear models. The R Journal 2(2), 20-28.

Shen, X., Alam, M., Fikse, F., and Rönnegård, L. (2013). A novel generalized ridge regression method for quantitative genetics. Genetics 193(4), 1255-1268.
}

\seealso{
\pkg{rrBLUP}\cr
\pkg{hglm}\cr
\pkg{glmnet}\cr
\pkg{spls}
}

\examples{
## Load example data ##
data("cubic")

## fit RRBLUP model ##
RRBLUP_model <- GSReModel(markers = Markers[1:200,],pheVal = pheData$DTT[1:200],modelMethods = "RRBLUP")

## fit SPLS model ##
SPLS_model <- GSReModel(markers = Markers[1:200,],pheVal = pheData$DTT[1:200],modelMethods = "SPLS")

## fit LASSO model ##
LASSO_model <- GSReModel(markers = Markers[1:200,],pheVal = pheData$DTT[1:200],modelMethods = "LASSO")

## fit bigRR model ##
bigRR_model <- GSReModel(markers = Markers[1:200,],pheVal = pheData$DTT[1:200],modelMethods = "bigRR")
}