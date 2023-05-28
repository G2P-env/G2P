\name{GSDataQC}
\alias{GSDataQC}

\title{
Data quality control of genomic selection
}
\description{
Function for summary and quality control of GS data.
}
\usage{
GSDataQC(markers, phenotype, impute = F, filter = F, NABound = 0.8,
         MAFBound = 0.05, imputeMethod = "median",hete = 1,silent =F)
}

\arguments{
  \item{markers}{(numeric, matrix) genotypic data, row represents sample and column represents SNP (feature). Missing (NA) alleles are allowed}
  \item{phenotype}{(data.frame) phenotypic data. The first column must be the sample name.}
  \item{impute}{(logical) if TRUE, imputation, default FALSE.}
  \item{filter}{(logical) if TRUE, filtering genotypic data with MAF and missing rate.}
  \item{NABound}{(numeric, [0,1]) max missing value percentage.}
  \item{MAFBound}{(numeric, [0,1]) min MAF percentage in each marker.}
  \item{imputeMethod}{(character) the method of imputation, "mean" or "median", default "median".}
  \item{hete}{(numeric) if genotypic matrix coded with numeric or transformed to numeric from letter, this parameter is used for MAF summary which define the numeric code of heterozygote, default 1 among 0,1,2.}
  \item{silent}{(logical) TRUE for hiding progress bar and FALSE for showing progress bar.}
}

\value{
a list of the data quality information.\cr
\item{Global}{overview of missing value in genotypic data.}
\item{resMarker}{genotypic data after processing.}
\item{MAF}{minor allele frequency (MAF) of each SNP.}
\item{NACountRow}{number of missing alleles from row perspective (sample).}
\item{NACountRowPercent}{missing rate from row perspective (sample).}
\item{NACountCol}{number of missing alleles from column perspective (feature SNP).}
\item{NACountColPercent}{missing rate from column perspective (feature SNP).}
\item{NAIndex}{index of all missing alleles.}
\item{phenotypeNAIdx}{index of missing phenotypes.}
}

\examples{
data("cubic")
## generate missing value
misIndex <- sample(1:length(Markers),1000000)
Markers[misIndex] <- NA

## GSDataQC, not impute ##
QCRes <- GSDataQC(markers = Markers, phenotype = pheData, impute = T)

## GSDataQC, impute with methods mean ##
QCResImpute <- GSDataQC(markers = Markers, phenotype = pheData, impute = T, 
                        imputeMethod = "mean")
}