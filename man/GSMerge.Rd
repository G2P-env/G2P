\name{GSMerge {G2P}}
\alias{GSMerge}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Results integration of GS methods
}
\description{
The function provides a strategy to integrate the results of two or more methods.
}
\usage{
GSMerge(predResMat, ratio, autoOptimize = F)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{predResMat}{(numeric, matrix) the prediction results of algorithms which you want to merge, the first column is the real value of trait.}
  \item{ratio}{(numeric,array) the weights of every algorithms.}
  \item{autoOptimize}{(logical) if true, two group results of methods from multiple-results will be selected and then get average (wight 1:1), default FALSE.}
}
\details{
The predResMat must including real value in first column, and if setting "autoOptimize = T", the count of algorithms must more than two.

In this function, if autoOptimize = T, the final two algorithms merge are selected from multi methods by following strategy:\cr
Firstly, compute the pearson's correlation of predResMat, choose the best correlation between real value and prediction scores, named \emph{method 1}.\cr
Secondly, choose the best correlation between method 1 and other methods, named \emph{method 2}.\cr
Finally, merge \emph{method 1} and \emph{method 2} with 1:1 (mean).

This auto-merge only based the evaluate metrics of \bold{pearson's correlation} evaluation
}
\value{
a matrix:\cr
involve observed value of trait, merge models and the merge result.
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
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
## load example data ##
data(cubic)

## Genotype-to-phenotype ##
G2Pres <- G2P(markers = Markers,data = pheData,trait = "EW",
              modelMethods = c("BayesA", "BRR", "RRBLUP",
                               "SPLS", "RFR"),trainIdx = 1:200,predIdx = 201:400)
                               
## merge 2 model prediction results ## 
GSMergeRes <- GSMerge(predResMat = G2Pres[,c("realPhenScore","RRBLUP","SPLS")], autoOptimize = F,ratio = c(3,7))

## merge 2 model prediction results automatically ## 
GSMergeRes <- GSMerge(predResMat = G2Pres, autoOptimize = T)
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
