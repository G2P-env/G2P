\name{GSEnsemble {G2P}}
\alias{GSEnsemble}
%- Also NEED an '\alias' for EACH other topic documented here.
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
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{predMat}{(numeric, matrix) the prediction results of algorithms which you want to merge, the first column is the real value of trait.}
  \item{nrandom}{(integer) the repeat number of stacking, default 10.}
  \item{evalMethods}{(character) ensemble base which evaluation methods.}
  \item{by}{(numeric,(0,1)) the radio window of ensemble, the smaller "by", the higher accuracy of ensemble. Default 0.1.}
  \item{evaluation}{(logical) if evaluate finalMat with evalMethods, default TRUE.}
  \item{topAlpha}{(numeric) the parameter of threshold evaluation methods, see also function evaluateGS. In this function, indicates the best ensemble base threshold when evalMethods is threshold methods.}
}
\details{
The function generates the integrative GS model by linearly combining multiple GS prediction results, the upgrade version of GSMerge. The weight of each individual GS model can be optimized by multiple iteration.
}
\value{
a list:\cr
\item{$BestWeight}{the best weight of methods in all repeat}
\item{$finalMat}{the final matrix cbind predMat with final ensemble score}
\item{$evalRes}{the evaluation results of finalMat with evalMethods}
\item{$weightMat}{a weight matrix including all repeats}
\item{$evalMat}{a evaluation results matrix including all repeats}
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\link{GSMerge}
\link{evaluateGS}
}
\examples{
## load example data ##
data(cubic)

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
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory (show via RShowDoc("KEYWORDS")):
% \keyword{ ~kwd1 }
% \keyword{ ~kwd2 }
% Use only one keyword per line.
% For non-standard keywords, use \concept instead of \keyword:
% \concept{ ~cpt1 }
% \concept{ ~cpt2 }
% Use only one concept per line.
