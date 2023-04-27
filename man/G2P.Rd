\name{G2P}
\alias{G2P}

\title{
Genotype to phenotype
}

\description{
The function fits genomic selection models, performs prediction and exports the prediction value of testing sets. 
}
\usage{
G2P(markers,data,fix = NULL,trait,trainIdx,predIdx,modelMethods ="BayesA",outputModel = FALSE,NAImpute = FALSE,
    nIter = 1500, burnIn = 500, thin = 5, saveAt = "", S0 = NULL, df0 =5, R2 = 0.5, weights = NULL,
    verbose = FALSE, rmExistingFiles = TRUE, groups=NULLimportance = FALSE,posPercentage = 0.4,
    BestIndividuals = c("top"),ntree = 500,nodesize = NULL,kernel = c("linear"),gamma = 1, cost = 2^(-9), 
    K = 8,eta = 0.7,select = "pls2",fit = "simpls",scale.x = FALSE,scale.y = FALSE,eps = 1e-4,maxstep = 100,
    parameters,alpha = 1,X = NULL,family = gaussian(link = identity), lambda = NULL,
    tol.err = 1e-6, tol.conv = 1e-8,epochs = 30, neurons = 4,
    ...)
}

\arguments{
  \item{markers}{(numeric, matrix) genootypic data, row represents sample and column represents SNP (feature).Genotypes generally be coded as {0,1,2};0 represents AA(homozygote),2 represents BB(homozygote) and 1 represents AB(heterozygote); missing (NA) alleles are not allowed.}
  \item{data}{(data.frame) phenotypic data and other information for data set.}
  \item{fix}{(character, array) names of data that used as fix effects, numeric, NAs are not permitted,}
  \item{testPheno}{(numeric) the phenotype value of test population individual, default NULL.}
  \item{trait}{(character) the trait name in data to fit model and perform genotype-to-phenotype prediction.}
  \item{trainIdx}{(numeric, array) index of training set in marker and data.}
  \item{predIdx}{(numeric, array)  index of testing/prediction set in marker and data.}
  \item{modelMethods}{(character) 16 alternative models to fit, including "BayesA", "BayesB", "BayesC", "BL", "BRR", "RKHS", "RRBLUP","LASSO", "SPLS", "SVC", "SVR", "RFR", "RFC", "RR", "RKHS", "BRNN"}
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
  \item{nodesize}{(integer) ramdonforest parameter Minimum size of terminal nodes. Setting this number larger causes smaller trees to be grown (and thus take less time). Note that the default values are different for classification (1) and regression (5).}
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
\details{
 \tabular{rlll}{
   \strong{Abbreviation}   \tab\tab  \strong{Full name} \cr
   \strong{Regression models} \cr
   BayesA   \tab\tab   Bayes A \cr
   BayesB   \tab\tab   Bayes B \cr
   BayesC   \tab\tab   Bayes C \cr
   BL   \tab\tab   Bayesion LASSO \cr
   BRR   \tab\tab   Bayesion ridge regression \cr
   RKHS  \tab\tab   Reproducing kernel Hilbert space \cr
   RR   \tab\tab   Ridge regression \cr
   RRBLUP   \tab\tab   Ridge regression best linear unbiased prediction \cr
   bigRR   \tab\tab   big ridge regression \cr
   LASSO   \tab\tab   Least absolute shrinkage and selection operator \cr
   SPLS   \tab\tab   Sparse partial least squares \cr
   SVR   \tab\tab   Support vector regression \cr
   RFR   \tab\tab   Random forest regression \cr
   BRNN  \tab\tab   Bayesian regularization neural networks \cr
   \strong{classification models} \cr
   SVC  \tab\tab   Support vector classification \cr
   RFC  \tab\tab   Random forest classification \cr
  }
}
\value{
a list including fitted model and prediction results (outputModel = TRUE).\cr
\item{model }{details of fitted model}
\item{predictRes }{a matrix of prediction results}

a matrix of prediction results (outputModel = FALSE).
}

\seealso{
\link[G2P]{GSmachine}\cr
\link[G2P]{GSReModel}\cr
\link[G2P]{fit.RR}\cr
\link[G2P]{fit.BRNN}\cr
\link[G2P]{fit.RKHS}\cr
\link[G2P]{predictGS}
}

\examples{
## load example data ##
data("cubic")

## perform genotype to phenotype prediction with all integrated methods ## 
C2Pres <- G2P(markers = Markers,data = pheData,trait = "DTT",
               modelMethods = c("BayesA", "BayesB", "BayesC", "BL", "BRR", 
                                "RKHS", "RRBLUP","LASSO", "SPLS", "SVC", "SVR", "RFR", "RFC",
                                "RR", "RKHS", "BRNN"),trainIdx = 1:200,predIdx = 201:400)
                                
## perform genotype to phenotype prediction with all integrated methods and output model ## 
C2Pres_outmodel <- G2P(markers = Markers,data = pheData,trait = "DTT",outputModel = T,
               modelMethods = c("BayesA", "BayesB", "BayesC", "BL", "BRR", 
                                "RKHS", "RRBLUP","LASSO", "SPLS", "SVC", "SVR", "RFR", "RFC",
                                "RR", "RKHS", "BRNN"),trainIdx = 1:200,predIdx = 201:400)
                                
## add other information as fix effects ## 
## fix effects not permit NA, so impute NA first.
data <- pheData
data <- apply(data[,-1], 2,function(x){x[which(is.na(x) == T)] <- mean(x,na.rm = T);x})
data <- as.data.frame(data)
data$EW <- pheData$EW

## perform genotype to phenotype prediction ## 
C2Pres_fix <- G2P(markers = Markers,data = data,trait = "EW",outputModel = T,fix = c("DTT","PH"),
                  modelMethods = c("BayesA", "BayesB", "BayesC", "BL", "BRR", 
                                   "RKHS", "RRBLUP","LASSO", "SPLS", "SVC", "SVR", "RFR", "RFC",
                                   "RR", "RKHS", "BRNN"),trainIdx = 1:200,predIdx = 201:400)
}

