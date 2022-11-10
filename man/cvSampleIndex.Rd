\name{cvSampleIndex {G2P}}
\alias{cvSampleIndex}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Generate sample index for training set and Testing Set
}
\description{
This function generates index of samples for N fold cross-validation.
}
\usage{
cvSampleIndex(sampleNum, cross = 5, seed = 1,randomSeed = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{sampleNum}{The number of samples needed to be partitioned into training and testing sets.}
  \item{cross}{The fold of cross validation.}
  \item{seed}{An integer used as the seed for data partition. The default value is 1.}
  \item{randomSeed}{Logical variable, default FALSE.}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
A list and each element including $trainIdx $testIdx and $cvIdx
\item{trainIdx}{The index of training samples.}
\item{testIdx}{The index of testing samples.}
\item{cvIdx}{The cross validation index.}
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
## leave-one out cross validation
a <- cvSampleIndex(sampleNum = nrow(Markers), cross = nrow(Markers), seed = 1)

## 5-fold cross validation
b <- cvSampleIndex(sampleNum = nrow(Markers), cross = 5, seed = 1)
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
