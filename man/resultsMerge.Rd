\name{resultsMerge {G2P}}
\alias{resultsMerge}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Transform prediction results list to matrix
}
\description{
The function provides a quick way to get prediction matix from prediction list.
}
\usage{
resultsMerge(predList)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{predList}{(list)  prediction results list of cross-validation.}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
a big matrix of prediction results.
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
## Load example data ##
data(GYSS)

## cross validation ##
predlist <- G2PCrossvalidation(cross = 10,seed = 1 , cpus = 3, markers = Markers,
               pheVal = phenotype, modelMethods = c("BayesA","BayesB","BayesC","rrBLUP", "RFC"),
               outputModel = FALSE)
resultMat <- resultsMerge(predlist)
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
