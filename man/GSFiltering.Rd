\name{GSFiltering}
\alias{GSFiltering}

\title{
Reduce genotypic data with Plink.
}
\description{
The function can raduce genotypic data by calling Plink.
}
\usage{
GSFiltering(input, output, inputformat, missingRate = 0.2, maf = 0.05, 
window = 100, step = 5, r2 = 0.1, workdir = "./", memory = 8)
}

\arguments{
  \item{input}{(character) name of original file, prefix for plink, full name for hapmap and vcf.}
  \item{output}{(character) name of original file, prefix for all.}
  \item{inputformat}{(character) file format of input file, "hmp" for hapmap,"vcf" VCF or "plink" for bed,bim,fam.}
  \item{missingRate}{(numeric [0-1]) filters out all variants with missing call rates exceeding the provided value to be removed}
  \item{maf}{(numeric [0-1]) filters out all variants with minor allele frequency below the provided threshold (default 0.05)}
  \item{window}{(numeric or character) window size in variant count (numeric) or kilobase (if the 'kb' modifier is present). }
  \item{step}{(numeric) variant count to shift the window at the end of each step}
  \item{r2}{(numeric) pairwise r2 threshold: at each step.}
  \item{workdir}{(character) set working directory}
  \item{memory}{(integer) set memory size, GB.}
}

\value{
several files after target filteration.\cr
\item{.hmp.txt }{genotypic data in hapmap format after filtering.}
\item{.prune.in }{SNP list of prune in.}
\item{.prune.out }{SNP list of prune out.}
\item{.log }{log file.}
}
\references{
https://www.cog-genomics.org/plink/1.9/filter#missing

https://www.cog-genomics.org/plink/1.9/ld#indep
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
GSFiltering(input = "test.hmp.txt", output = "filter_res", 
            inputformat = "hmp", missingRate = 0, 
            maf = 0.05,window = 100,step = 5,r2 = 0.1)
}