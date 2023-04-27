\name{GSRead}
\alias{GSRead}

\title{
Fast read of the larger genotypic data
}
\description{
Fast read of the larger genotypic data in PLINK HapMap and VCF format, and convert it to numeric format (0 1 2) by calling Plink.
}
\usage{
GSRead(input, output, inputformat, workdir = "./", memory = 8)
}

\arguments{
  \item{input}{(character) name of original file, prefix for plink, full name for hapmap and vcf.}
  \item{output}{(character) name of original file, prefix for all.}
  \item{inputformat}{(character) file format of input file, "hmp" for hapmap,"vcf" VCF or "plink" for bed,bim,fam.}
  \item{workdir}{(character) set working directory}
  \item{memory}{(integer) set memory size, GB.}
}

\value{
several files after target filteration.\cr
\item{*_temp* }{genotypic data in other formats when format transformation is done.}
\item{.raw }{recode file in numeric (0,1,2) format by Plink.}
\item{.log }{log file.}
}
\references{
https://www.cog-genomics.org/plink/1.9/data#recode
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
write.table(hmp,file = "test.hmp.txt",row.names = F,col.names = T,sep = "\t",quote = F)

## filtering 
G <- GSRead(input = "test.hmp.txt", output = "read_tmp", 
     inputformat = "hmp", workdir = "./", memory = 8)

## show genotypic in numeric
dim(G)
G[1:10,1:10]
}