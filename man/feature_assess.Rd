\name{feature_assess {G2P}}
\alias{feature_assess}
\title{
Feature assessment and selection
}

\description{
This function scores each marker,so that reduce the data dimension and perform feature selection.(Methods including Gini,Accuracy and rrBLUP).
}
\usage{
feature_assess(markers, phenotype, method = c("RRBLUP", "Gini", "Accuracy"), 
               ntree = 500, importance = TRUE, posPercentage = 0.40, 
               BestIndividuals = c("top", "middle", "buttom"))
}

\arguments{
  \item{markers}{(numeric, matrix) row represents sample and column represent SNP information (feature).Genotypes should be coded as {0,1,2};0 represents AA(homozygote),2 represents BB(homozygote) and 1 represents AB(heterozygote);missing alleles (NAs) are not allowed.}
  \item{phenotype}{(numeric) a numeric arrary, phenotype value of each individual.}
  \item{method}{(character) the method of feature selction, including "Gini" "Accuracy" "rrBLUP", default "rrBLUP".}
  \item{ntree}{(numeric) the number of random forest decision tree, default 500.}  
  \item{importance}{(logical) whether the results of variable importance,default TRUE.}
  \item{posPercentage}{(numeric,[0,1]) phenotype of extreme individuals which expected, default 0.4.}
  \item{BestIndividuals}{(character) the position of expected phenotypes in whole phenotypic dataset."top","buttom" or "middle",defult "top".}
}

\value{
A numeric arrary, regarding to score of each position of SNPs
}
\references{
Endelman, J. B. (2011). Ridge regression and other kernels for genomic selection with R package rrBLUP. The plant genome, 4(3).

Breiman, L. (2001), Random Forests, Machine Learning 45(1), 5-32.
}

\examples{
data("cubic")

Gini_selection <- feature_assess(markers = Markers[1:100,], phenotype = pheData$DTT[1:100], method = "Gini", 
                                 ntree = 500, importance = TRUE, posPercentage = 0.40,
                                 BestIndividuals = "top")

## feature selection with Acc ##
Acc_selection <- feature_assess(markers = Markers[1:100,], phenotype = pheData$DTT[1:100], method = "Accuracy", 
                                ntree = 500, importance = TRUE, posPercentage = 0.40,
                                BestIndividuals = "top")

## feature selection with rrBLUP ##
rrBLUP_selection <- feature_assess(markers = Markers[1:100,], phenotype = pheData$DTT[1:100], method = "RRBLUP", 
                                   posPercentage = 0.40, BestIndividuals = "top")
}