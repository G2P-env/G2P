\name{GSTransForm}
\alias{GSTransForm}

\title{
Genotypic data format conversion
}
\description{
By calling other external software like VCFtools Plink TASSEL etc for the transition of genotype format.
}
\usage{
GSTransForm(pattern, input, output, workdir = "./", memory = 8 )
}

\arguments{
  \item{pattern}{(character) pattern of transformation, including "hmp2plink","plink2vcf","vcf2hmp","vcf2plink","plink2hmp", and "hmp2vcf".}
  \item{input}{(character) name of original file, prefix for plink, full name for hapmap and vcf.}
  \item{output}{(character) name of original file, prefix for all.}
  \item{workdir}{(character) set working directory}
  \item{memory}{(integer) set memory size }
}

\value{
one to several files in target format.
}

\note{
The function only be accessible for G2P container.
}

\seealso{
https://www.cog-genomics.org/plink/

https://www.maizegenetics.net/tassel

https://vcftools.sourceforge.net/
}

\examples{
## load example data ##
library(G2P)
data("cubic_hmp")
setwd("./")
write.table(hmp,file = "test.hmp.txt",row.names = F,col.names = T,sep = "\t",quote = F)

## transformation ##
GSTransForm(pattern="hmp2plink",input="test.hmp.txt",output="h2p",workdir="./")
}