\name{sampleClassify}
\alias{sampleClassify}

\title{
Generate positive and negative samples for classification methods.
}

\description{
The function can be used to generate positive and negative samples for training.The positive samples represent the excellent individuals which's breeding values we expect to obtain in research.And the negative samples represent the lower breeding values of individuals.
}
\usage{
sampleClassify(phenotype, posPercentage = 0.4,
               BestIndividuals = c("top", "middle", "buttom"))
}

\arguments{
  \item{phenotype}{(numeric) breeding values or phenotypic values of individuals.}
  \item{posPercentage}{(numeric,[0,1]) phenotype of extreme individuals which expected, default 0.4.}
  \item{BestIndividuals}{(character)the position of expected phenotype in whole phenotypic data set."top","buttom" or "middle",default "top".}
}

\value{
A list of row number regarding to positive and negative samples.
\item{posSampleIndex }{index of positive samples.}
\item{negSampleIndex }{index of negative samples.}

}

\examples{
data("cubic")

## percentage of positive samples is 0.4 ##
sampleCly <- sampleClassify(pheData$DTT, posPercentage = 0.4, BestIndividuals = "top")
}