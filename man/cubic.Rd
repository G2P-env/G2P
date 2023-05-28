\name{cubic}
\alias{cubic}
\docType{data}

\title{
Example data for G2P
}
\description{
Maize CUBIC population, after filterion (-indep-pairwise 100 5 0.1 -indep-pairwise 100 5 0.1), remaining 1,404 lines and 9,286 high quality SNPs, coded by 0,1,2. \cr

Phenotypic data involves tree traits "DTT" for days to tasseling , "PH" for plant height and "EW" for ear weight.
}
\usage{
  data("cubic")
}
\format{
  \describe{
    \item{\code{Markers}}{a numeric matrix of genotypes.}
    \item{\code{pheData}}{a dataframe of tree traits phenotypes.}
  }
}

\source{
http://cubicmaize.hzau.edu.cn/
}

\references{
Liu, HJ., Wang, X., Xiao, Y. et al. CUBIC: an atlas of genetic architecture promises directed maize improvement. Genome Biol 21, 20 (2020). https://doi.org/10.1186/s13059-020-1930-x

Jingyun Luo, Chengcheng Wei, Haijun Liu, Shikun Cheng, Yingjie Xiao, Xiaqing Wang, Jianbing Yan, Jianxiao Liu, MaizeCUBIC: a comprehensive variation database for a maize synthetic population, Database, Volume 2020, 2020, baaa044, https://doi.org/10.1093/database/baaa044
}

\examples{
data("cubic")
}