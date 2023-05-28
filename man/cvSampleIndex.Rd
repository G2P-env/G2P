\name{cvSampleIndex}
\alias{cvSampleIndex}

\title{
Generate sample index for training set and test set
}

\description{
This function generates index of samples for N fold cross-validation.
}
\usage{
cvSampleIndex(sampleNum, cross = 5, seed = 1,randomSeed = FALSE)
}
\arguments{
  \item{sampleNum}{(integer) the number of samples needed to be partitioned into training and test set.}
  \item{cross}{(integer) the fold of cross validation.}
  \item{seed}{(integer) the seed for data partition. The default value is 1.}
  \item{randomSeed}{(logical) variable, default FALSE.}
}

\value{
A list including "trainIdx" "testIdx" and "cvIdx"
\item{trainIdx}{index of training samples.}
\item{testIdx}{index of test samples.}
\item{cvIdx}{cross validation index.}
}

\examples{
## Load example data ##
data("cubic")
## leave-one out cross validation
a <- cvSampleIndex(sampleNum = nrow(Markers), cross = nrow(Markers), seed = 1)

## 5-fold cross validation
b <- cvSampleIndex(sampleNum = nrow(Markers), cross = 5, seed = 1)
}
