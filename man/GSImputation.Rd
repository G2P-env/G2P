\name{GSImputation}
\alias{GSImputation}

\title{
Imputation of genotypic data with Beagle.
}
\description{
The function can execute imputation of genotypic data by calling Beagle.
}
\usage{
GSImputation(input, output, inputformat, workdir = "./", 
             memory = 8,threads = 1,window = 40,overlap = 2)
}

\arguments{
  \item{input}{(character) name of original file, prefix for plink, full name for hapmap and vcf.}
  \item{output}{(character) name of original file, prefix for all.}
  \item{inputformat}{(character) file format of input file, "hmp" for hapmap,"vcf" VCF or "plink" for bed,bim,fam.}
  \item{memory}{(integer) set memory size, GB}
  \item{threads}{(integer) specify the number of computational threads that will be used.}
  \item{window}{(numeric) specify the cM length of each sliding window (default: 40).}
  \item{overlap}{(numeric) specify the cM length of overlap between adjacent sliding windows(default: 2).}
  \item{workdir}{(character) set working directory}
}

\value{
several files after target filteration.\cr
\item{.vcf }{genotypic data after imputation in VCF format.}
\item{.log }{log file.}
}
\references{
http://faculty.washington.edu/browning/beagle/beagle.html

B L Browning, X Tian, Y Zhou, and S R Browning (2021) Fast two-stage phasing of large-scale sequence data. Am J Hum Genet 108(10):1880-1890. doi:10.1016/j.ajhg.2021.08.005

B L Browning, Y Zhou, and S R Browning (2018). A one-penny imputed genome from next generation reference panels. Am J Hum Genet 103(3):338-348. doi:10.1016/j.ajhg.2018.07.015 
}

\note{
The function only be accessible for G2P container.
}

\seealso{
\link[G2P]{GSTransForm}
}
\examples{
## load example data ##
library(G2P)
data("cubic_hmp")
setwd("./")

## generating NA to data
hmp_NA <- as.matrix(hmp[,-c(1:11)])
hmp_NA[sample(1:length(hmp_NA),1000000)] <- "NN"
hmp_NA <- cbind(hmp[,1:11],hmp_NA)
write.table(hmp_NA,file = "test_NA.hmp.txt",row.names = F,col.names = T,sep = "\t",quote = F)

## imputation 
GSImputation(input = "test_NA.hmp.txt", output = "impute_res_NA",
             workdir = "./", inputformat = "hmp")
}