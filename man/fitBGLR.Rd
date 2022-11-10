\name{fit.BGLR {G2P}}
\alias{fit.BGLR}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Fit Bayesian generalized linear regression
}
\description{
The function fits various types of parametric and semi-parametric Bayesian regressions.
}
\usage{
fit.BGLR(trainMarkerMat, trainPheVal, predictMarkerMat, modelMethods,
         trainX = NULL, predX = NULL, outputModel = FALSE, nIter = 1500, 
         burnIn = 500, thin = 5,verbose = FALSE,
         saveAt = "", S0 = NULL, df0 =5, R2 = 0.5, weights = NULL,
         rmExistingFiles = TRUE, groups=NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{trainMarkerMat}{(numeric, matrix) genotypic data of training set. row represents sample and column represents SNP (feature).Genotypes generally be coded as {0,1,2};0 represents AA(homozygote),2 represents BB(homozygote) and 1 represents AB(heterozygote); missing (NA) alleles are not allowed.}
  \item{trainPheVal}{(numeric, array) the phenotypic value of each sample.}
  \item{predictMarkerMat}{(numeric, matrix) genotypic data of testing set. row represents sample and column represents SNP (feature).Genotypes generally be coded as {0,1,2};0 represents AA(homozygote),2 represents BB(homozygote) and 1 represents AB(heterozygote); missing (NA) alleles are not allowed.}
  \item{modelMethods}{(character) the model to fit. "BayesA", "BayesB", "BayesC", "BL", "BRR".}
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
  \item{outputModel}{(logical) if TRUE, return the list of trained model and prediction results, default FALSE.}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
a list including fitted model and prediction results (outputModel = TRUE).\cr
\item{$model }{details of fitted model}
\item{$predictRes }{prediction results of testing set}

a array of prediction results (outputModel = FALSE).
%% ...
}
\references{
Pérez, P., & de los Campos, G. (2014). BGLR: a statistical package for whole genome regression and prediction. Genetics, 198(2), 483-495.
}
\author{
%%  ~~who you are~~
}
\note{
For more details, see http://genomics.cimmyt.org/BGLR-extdoc.pdf
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\pkg{BGLR}
}
\examples{
## Load example data ##
data(cubic)

## Fit Bayes A model and output model details ##
BayesA_model <- fit.BGLR(trainMarkerMat = Markers[1:200,], predictMarkerMat = Markers[201:400,],
                         trainPheVal = pheData$DTT[1:200], modelMethods = "BayesA", outputModel = TRUE)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory (show via RShowDoc("KEYWORDS")):
% \keyword{ ~kwd1 }
% \keyword{ ~kwd2 }
% Use only one keyword per line.
% For non-standard keywords, use \concept instead of \keyword:
% \concept{ ~cpt1 }
% \concept{ ~cpt2 }
% Use only one concept per line.
