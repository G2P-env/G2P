\name{G2PEvaluation {G2P}}
\alias{G2PEvaluation}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Models evaluation
}
\description{
The function consolidates 13 evaluation metrics to evaluete the performance of genomic selection models.
}
\usage{
G2PEvaluation(realScores, predScores, Probability = FALSE, evalMethod = "RE",
           Beta = 1, BestIndividuals = "top", topAlpha = 1:90,probIndex = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{realScores}{(numeric,vector) observation values of validation data.}
  \item{predScores}{(numeric,vector or matrix) prediction values from genomic selection models.}
  \item{evalMethod}{(character) alternative evaluation metrics, including "pearson", "kendall", "spearman", "MSE", "R2", "RE", "Kappa", "AUC", "AUCpr", "accuracy", "F1", "meanNDCG", "NDCG".}
  \item{Beta}{(numeric) the parameter of F-measure or F-score, when beta = 1, meaning F1, Default F1.}
  \item{BestIndividuals}{(character) the position of expected phenotype in whole phenotypic data set."top","buttom" or "middle",default "top".}
  \item{topAlpha}{(numeric, array) the proportion of considering excellent samples, default 1:90.}
  \item{globalAlpha}{(logical) if global metrics (pearson, kendall, spearman, MSE and R2) evaluate by threshold (topAlpha),default FALSE.}
  \item{probIndex}{(integer) indicate the column index which prediction result is probability. For example, classification model "RFC" locates the fifth column in prediction matrix and then probIndex = 5.}
}
\details{
 \tabular{rlll}{
   \strong{Abbreviation}   \tab\tab  \strong{Full name} \cr
   pearson   \tab\tab   Pearson’s correlation coefficient \cr
   kendall   \tab\tab   Kendall’s rank correlation coefficient \cr
   spearman   \tab\tab   Spearman's rank correlation \cr
   MSE   \tab\tab   Mean squared error \cr
   R2   \tab\tab   Coefficient of determination \cr
   RE   \tab\tab   Relative efficiency \cr
   Accuracy   \tab\tab   Accuracy \cr
   Kappa   \tab\tab   Cohen's kappa coefficient \cr
   AUC   \tab\tab   Area under the receiver operating characteristic curve \cr
   AUCpr   \tab\tab   Area under the precision-recall curve \cr
   F1   \tab\tab   F-score(beta = 1) \cr
   meanNDCG   \tab\tab   Mean normalized discounted cumulative gain \cr
   NDCG   \tab\tab   Normalized discounted cumulative gain \cr
  }
}
\value{
a list inculding evaluation results with methods which user selected.
\item{corMethods }{a matrix for global metrics}

    \item{RE }{a matrix for threshold metric RE.}
    \item{Kappa }{a matrix for threshold metric Kappa.}
    \item{... }{a matrix for other threshold metrics user selected.}
}
\references{
Blondel M, Onogi A, Iwata H, Ueda N. 2015. A ranking approach to genomic selection. PloS One, 10(6): e0128570
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
data(cubic)
########## prediction ############
C2Pres <- G2P(markers = Markers,data = pheData,trait = "EW",
              modelMethods = c("BayesA", "BRR", "RRBLUP",
                               "SPLS", "SVC", "RFR"),trainIdx = 1:200,predIdx = 201:400)
########## evaluation ############
evalres <- G2PEvaluation(realScores = C2Pres[,1], predScores = C2Pres[,2:7], 
                           evalMethod = c("pearson", "kendall","spearman","RE","Kappa",
                                          "AUC","AUCpr","NDCG","meanNDCG",
                                          "MSE","R2","F1","accuracy"),topAlpha = 1:90, probIndex = 5)
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
