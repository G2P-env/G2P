\name{resultsMerge}
\alias{resultsMerge}

\title{
Transform prediction results list to matrix
}
\description{
The function provides a quick way to get prediction matix from prediction list.
}
\usage{
resultsMerge(predList)
}

\arguments{
  \item{predList}{(list)  prediction results list of cross-validation.}
}

\value{
a big matrix of prediction results.
}

\examples{
## Load example data ##
data("cubic")

## cross validation ##
predlist <- G2PCrossValidation(cross = 5,seed = 1, cpus = 1, times = 2, markers = Markers[1:200,],
                data = pheData[1:200,], trait = "DTT" modelMethods = c("rrBLUP", "spls"),
                outputModel = FALSE)
resultMat <- resultsMerge(predlist)
}