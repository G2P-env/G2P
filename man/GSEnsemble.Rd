\name{GSEnsemble}
\alias{GSEnsemble}

\title{
Results integration of multiple GS methods
}

\description{
The function provides a strategy to integrate the results of multiple GS methods.
}
\usage{
GSEnsemble(predMat, nrandom = 10, evalMethods, by = 0.1, evaluation = T, topAlpha = 15, 
           ...)
}

\arguments{
  \item{predMat}{(numeric, matrix) prediction results of algorithms which you want to merge, the first column is the real value of trait.}
  \item{nrandom}{(integer) repeat number of stacking, default 10.}
  \item{evalMethods}{(character) ensemble base which evaluation methods.}
  \item{by}{(numeric,(0,1)) window of ensemble, the smaller "by", the higher accuracy of ensemble. Default 0.1.}
  \item{evaluation}{(logical) if evaluate finalMat with evalMethods, default TRUE.}
  \item{topAlpha}{(numeric) parameter of threshold evaluation methods, see also function evaluateGS. In this function, indicates the best ensemble base threshold when evalMethods is threshold methods.}
}
\details{
The function generates the integrative GS model by linearly combining multiple GS prediction results, the upgrade version of GSMerge. The weight of each individual GS model can be optimized by multiple iteration.
}
\value{
a list:\cr
\item{BestWeight}{best weight of methods in all repeat}
\item{finalMat}{final matrix cbind predMat with final ensemble score}
\item{evalRes}{evaluation results of finalMat with evalMethods}
\item{weightMat}{a weight matrix including all repeats}
\item{evalMat}{a evaluation results matrix including all repeats}
}

\seealso{
\link{GSMerge}
\link{evaluateGS}
}
\examples{
## load example data ##
data("cubic")

## Genotype-to-phenotype ##
G2Pres <- G2P(markers = Markers,data = pheData,trait = "EW",
              modelMethods = c("BayesA", "BRR", "RRBLUP",
                               "SPLS", "RFR"),trainIdx = 1:200,predIdx = 201:400)

## ensemble more than 2 prediction results of models automatically ## 
GSEnsembleRes <- GSEnsemble(predMat = G2Pres, nrandom = 10, evalMethods = "pearson",
                          by = 0.1, evaluation = T)
## new prediction of genotype-to-phenotype ##
G2PresNew <- G2P(markers = Markers,data = pheData,trait = "EW",
              modelMethods = c("BayesA", "BRR", "RRBLUP",
                               "SPLS", "RFR"),trainIdx = 1:200,predIdx = 401:600)

## apply GSEnsemble rotio to new prediction ##                       
GSEnsembleModelRes <- GSMerge(predResMat = G2PresNew,ratio = GSEnsembleRes$BestWeight, autoOptimize = F)                      
}

