\name{cvSampleIndex}
\alias{cvSampleIndex}

\title{
Generate sample index for training set and Testing Set
}

\description{
This function generates index of samples for N fold cross-validation.
}
\usage{
cvSampleIndex(sampleNum, cross = 5, seed = 1,randomSeed = FALSE)
}
\arguments{
  \item{sampleNum}{The number of samples needed to be partitioned into training and testing sets.}
  \item{cross}{The fold of cross validation.}
  \item{seed}{An integer used as the seed for data partition. The default value is 1.}
  \item{randomSeed}{Logical variable, default FALSE.}
}

\value{
A list and each element including "trainIdx" "testIdx" and "cvIdx"
\item{trainIdx}{The index of training samples.}
\item{testIdx}{The index of testing samples.}
\item{cvIdx}{The cross validation index.}
}

\examples{
## Load example data ##
data("cubic")
## leave-one out cross validation
a <- cvSampleIndex(sampleNum = nrow(Markers), cross = nrow(Markers), seed = 1)

## 5-fold cross validation
b <- cvSampleIndex(sampleNum = nrow(Markers), cross = 5, seed = 1)
}
