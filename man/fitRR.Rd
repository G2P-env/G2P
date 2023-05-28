\name{fit.RR}
\alias{fit.RR}

\title{
Fitting ridge regression model and performing prediction
}
\description{
The function fits ridge regression model and export the prediction value of test set. 
}

\usage{
fit.RR(trainMarkerMat, trainPheVal, predictMarkerMat, NAImpute = F,
       trainX = NULL,predX = NULL,cpus = 1,outputModel = FALSE)
}

\arguments{
  \item{trainMarkerMat}{(numeric, matrix) genotypic data of training set. row represents sample and column represents SNP (feature).Genotypes generally be coded as {0,1,2};0 represents AA(homozygote),2 represents BB(homozygote) and 1 represents AB(heterozygote); missing (NA) alleles are not allowed.}
  \item{trainPheVal}{(numeric) the phenotypic value of each sample.}
  \item{predictMarkerMat}{(numeric, matrix) genotypic data of test set. row represents sample and column represents SNP (feature).Genotypes generally be coded as {0,1,2};0 represents AA(homozygote),2 represents BB(homozygote) and 1 represents AB(heterozygote), the coding pattern must be consistent with training data. missing (NA) alleles are not allowed.}
  \item{NAImpute}{(logical) if true, the missing values in phenotype will be imputed by mean otherwise be removed.}
  \item{trainX}{additional data inputs of training sets as fixed effects instead of SNPs inputs only.}
  \item{predX}{additional data inputs of test set as fixed effects instead of SNPs inputs only.}
  \item{cpus}{(integer) number of cpu cores to use for calculations (only available in UNIX-like operating systems), default 1.}
  \item{outputModel}{(logical) if true, return the list of training model.}
}

\value{
a list including fitted model and prediction results (outputModel = TRUE).\cr
\item{model }{details of fitted model}
\item{predictRes }{prediction results of test set}

an array of prediction results (outputModel = FALSE).
}
\references{
Endelman, J.B. 2011. Ridge regression and other kernels for genomic selection with R package rrBLUP. Plant Genome 4:250-255. doi: 10.3835/plantgenome2011.08.0024
}

\seealso{
\pkg{rrBLUP}
}
\examples{
## Load example data ##
data("cubic")

## use RR model to modeling and predict ##
RR_res <- fit.RR(trainMarkerMat = Markers[1:200,], trainPheVal = pheData$DTT[1:200],
                 predictMarkerMat = Markers[201:400,],outputModel = T)
}