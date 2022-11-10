\name{sampleClassify {G2P}}
\alias{sampleClassify}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Generate positive and negative samples for training
}
\description{
The function can be used to generate positive and negative samples for training.The positive samples represent the excellent individuals which's breeding values we expect to obtain in research.And the negative samples represent the lower breeding values of individuals.
}
\usage{
sampleClassify(phenotype, posPercentage = 0.4,
               BestIndividuals = c("top", "middle", "buttom"))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{phenotype}{(numeric) breeding values or phenotypic values of individuals.}
  \item{posPercentage}{(numeric,[0,1]) phenotype of extreme individuals which expected, default 0.4.}
  \item{BestIndividuals}{(character)the position of expected phenotype in whole phenotypic data set."top","buttom" or "middle",default "top".}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
A list of row number regarding to positive and negative samples.
\item{posSampleIndex }{Index of positive samples.}
\item{negSampleIndex }{Index of negative samples.}
%% ...
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
data(cubic)

## percentage of positive samples is 0.4 ##
sampleCly <- sampleClassify(phenotype, posPercentage = 0.4, BestIndividuals = "top")
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
