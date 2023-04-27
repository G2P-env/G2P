\name{transLetter2number}
\alias{transLetter2number}

\title{
Transformate genotypic data from letter to number.
}
\description{
The function was designed for transforming genotypic data from letter to number. The input just a dataframe or matrix with SNP markers in columns and individuals in rows. Only genotypic data, NO extra columns of plant names or extra rows of markers information.
}
\usage{
transLetter2number(G,maf=0,mr=0,silent = FALSE,ref.alt = NULL,
                   impute = TRUE,imputeMethod = "median")
}

\arguments{
\item{G}{(character, matrix or dataframe) a matirx or dataframe with markers in columns and individuals in rows, NAs are allowed.}
\item{maf}{(numeric, [0,1]) minor allele frequency used to filter the SNP markers, the default is zero.}
\item{mr}{(numeric, [0,1]) missing rate of allele used to filter the SNP markers, the default is zero.}
\item{silent}{if a progress bar should be drawn for each step of the conversion. The default is FALSE.}
\item{impute}{(logical) if TRUE, imputation, default FALSE.}
\item{imputeMethod}{(character) the method of imputation, "mean" or "median", default "median".}
\item{ref.alt}{(character, matrix or dataframe) a matrix or dataframe with the information of reference alleles to be used for the transformation. The matrix or dataframe should have as many rows as markers with reference alleles and with 2 columns, being the first column of reference allele (ref) and the second column of alt allele (alt). If this information not be provided,the function will refer the MAF results.}
}

\value{
a list:
\item{G.n}{numeric matrix with markers in columns and individuals in rows. Coding by 0:ref homozygote, 1:heterozygote, 2:alt homozygote.} 
\item{maf.info}{the information of ref, alt, AF and MAF.\cr}
}

\seealso{
\link[G2P]{GSDataQC}
}

\examples{
## load hmp data
data("cubic_hmp")

## remove annotation infomation and inverse
G.hmp <- hmp
rownames(G.hmp) <- G.hmp[,1]
G.hmp <- G.hmp[,-c(1:11)]
G.hmp <- t(G.hmp)

## transform letter to numeric
G.n <- transLetter2number(G.hmp,maf=0,mr=0,silent = FALSE,ref.alt = NULL)
}