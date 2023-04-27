\name{G2PCrossValidation}
\alias{G2PCrossValidation}

\title{
Cross-validation for genomic selection
}
\description{
The function achieve cross-validation for integrated GS models.
}
\usage{
G2PCrossValidation(cvSampleList = NULL,cross = 10,times = 1,seed = 1,cpus = 1,markers,data,trait, 
                   modelMethods ="SVC", outputModel = FALSE,NAImpute = FALSE
                   nIter = 1500, burnIn = 500, thin = 5, 
                   saveAt = "", S0 = NULL, df0 =5, R2 = 0.5, weights = NULL,
                   verbose = FALSE, rmExistingFiles = TRUE, groups=NULL, importance = FALSE,
                   posPercentage = 0.4, BestIndividuals = c("top"), ntree = 500, nodesize = NULL, kernel = c("linear"), gamma = 1, cost = 2^(-9),
                   K = 8, eta = 0.7, select = "pls2", fit = "simpls", scale.x = FALSE, scale.y = FALSE, eps = 1e-4, trace = FALSE, maxstep = 100,
                   alpha = 1,family = gaussian(link = identity), lambda = NULL, tol.err = 1e-6, tol.conv = 1e-8,
                   epochs = 30, neurons = 4,
                   ...)
}

\arguments{
  \item{cvSampleList}{(list) the specified cross-validation combination, or a list from function \code{\link{cvSampleIndex}}.Default NULL and execute a random grouping.}
  \item{cross}{(integer) the fold of N-fold cross-validation. For example cross = 10 for 10-fold cross-validation.}
  \item{times}{(integer) the times of repetition. For example, cross = 10 and times = 5 for five times 10-fold cross-validation.}
  \item{seed}{(integer) random seed for sampling of cross-validation.}
  \item{cpus}{(integer) number of cpu cores to parallelizing cross-validation(only available in UNIX-like operating systems), default 1.}
  \item{markers}{(numeric, matrix) row represents sample and column represents SNP (feature).Genotypes generally be coded as {0,1,2};0 represents AA(homozygote),2 represents BB(homozygote) and 1 represents AB(heterozygote); missing (NA) alleles are not allowed.}
  \item{data}{(data.frame) phenotypic data and other information for data set.}
  \item{trait}{(character) name of aim trait to perform cross-validation.}
  \item{modelMethods}{(character) 16 alternative models to fit, including "BayesA, "BayesB, "BayesC", "BL", "BRR", "RKHS", "RRBLUP","LASSO", "SPLS", "SVC", "SVR", "RFR", "RFC", "RR", "RKHS", "BRNN"}
  \item{outputModel}{(logical) if TRUE, return the list of trained model and prediction results, default FALSE.}
  \item{NAImpute}{(logical) if true, the missing values in phenotype will be imputed by mean otherwise be removed.}
  \item{nIter,burnIn,thin}{(integer) the number of iterations, burn-in and thinning,default nIter 7000,burnIn 500,thin 5.}
  \item{saveAt}{(string) this may include a path and a pre-fix that will be added to the name of the files that are saved as the program runs,default "".}  
  \item{S0,df0}{(numeric) the scale parameter for the scaled inverse-chi squared prior assigned to the residual variance, only used with Gaussian outcomes.Default S0 NULL,df0 = 5.}
  \item{R2}{(numeric, (0,1)) the proportion of variance that one expects, a priori, to be explained by the regression. Only used if the hyper-parameters are not specified.Defult 0.5}
  \item{weights}{(numeric) a vector of weights, may be NULL. If weights is not NULL, the residual variance of each data-point is set to be proportional to the square of the weight. Only used with Gaussian outcomes.}
  \item{verbose}{(logical) if TRUE the iteration history is printed, default FALSE.}
  \item{rmExistingFiles}{(logical) if TRUE, removes existing output files from previous runs, default TRUE.}
  \item{groups}{(factor) a vector of the same length of y that associates observations with groups, each group will have an associated variance component for the error term.}
  \item{ntree}{(integer) randomforest parameter. Number of trees to grow. This should not be set to too small a number, to ensure that every input row gets predicted at least a few times,default 500.}
  \item{nodesize}{(integer) randomforest parameter Minimum size of terminal nodes. Setting this number larger causes smaller trees to be grown (and thus take less time). Note that the default values are different for classification (1) and regression (5).}
  \item{posPercentage}{(numeric,[0,1]) phenotype of extreme individuals which expected, default 0.4.}
  \item{BestIndividuals}{(character) the position of expected phenotype in whole phenotypic data set."top","buttom" or "middle",default "top".}
  \item{kernel}{(numeric) svm parameter the kernel used in training and predicting. You might consider changing some of the following parameters, depending on the kernel type.(linear,polynomial,sigmoid,radial). Default "linear".}
  \item{gamma}{(numeric)svm parameter parameter needed for all kernels except linear, default 1.}
  \item{cost}{(numeric)svm parameter cost of constraints violation, default: 2^(-9), it is the 'C'-constant of the regularization term in the Lagrange formulation.}
  \item{K}{(integer) SPLS parameter: number of hidden components.}
  \item{eta}{(numeric) SPLS parameter: thresholding parameter. eta should be between 0 and 1.}
  \item{select}{(character) SPLS parameter: PLS algorithm for variable selection. Alternatives are "pls2" or "simpls". Default is "pls2".} 
  \item{fit}{(character) SPLS parameter: PLS algorithm for model fitting. Alternatives are "kernelpls", "widekernelpls", "simpls", or "oscorespls". Default is "simpls".}
  \item{scale.x}{(character) SPLS parameter: scale predictors by dividing each predictor variable by its sample standard deviation?}
  \item{scale.y}{(character) SPLS parameter: scale responses by dividing each response variable by its sample standard deviation?}
  \item{eps}{(character) SPLS parameter: an effective zero. Default is 1e-4.}
  \item{maxstep}{(integer) SPLS parameter: maximum number of iterations when fitting direction vectors. Default is 100.}
  \item{alpha}{(numeric) LASSO model parameter: the elasticnet mixing parameter.Detail in glmnet.}
  \item{family}{the distribution family of y, see help('family') for more details. }
  \item{lambda}{the shrinkage parameter determines the amount of shrinkage. Default is NULL meaning that it is to be estimated along with other model parameters.}
  \item{tol.err}{internal tolerance level for extremely small values; default value is 1e-6.}
  \item{tol.conv}{tolerance level in convergence; default value is 1e-8.}
  \item{neurons}{(integer) indicates the number of neurons,defult 4.}
  \item{epochs}{(integer) maximum number of epochs(iterations) to train, default 30.}
}

\value{
a list:
The prediction results of input GS methods with cross validation.
\item{CV1 }{prediction matrix of first time CV}
\item{CV2 }{prediction matrix of second time CV}
\item{CV\emph{n}... }{prediction matrix of \emph{n}th time CV}
}

\seealso{
\code{\link{cvSampleIndex}}\cr
\code{\link{G2P}}\cr
}
\examples{
## load example data ##
data("cubic")

## perform a 2 times 5-fold cross-validation  ##
predlist <- G2PCrossValidation(cross = 5,seed = 1, cpus = 1, times = 2, markers = Markers[1:200,],
                data = pheData[1:200,], trait = "DTT" modelMethods = c("rrBLUP", "spls"),
                outputModel = FALSE)
}