\name{cubic_hmp}
\alias{cubic_hmp}
\docType{data}

\title{
    Example data in Hapmap format for G2P
}

\description{
Genotypic data of maize CUBIC population (1428 lines with 14,349,237 SNP), after filterion (-indep-pairwise 100 5 0.1 -indep-pairwise 100 5 0.1), remaining 1,428 lines (24 founder + 1,404 progeny) and 9,286 high quality SNPs in standard Hapmap format. \cr
}

\usage{
    data("cubic_hmp")
}

\source{
http://cubicmaize.hzau.edu.cn/
}

\references{
Liu, HJ., Wang, X., Xiao, Y. et al. CUBIC: an atlas of genetic architecture promises directed maize improvement. Genome Biol 21, 20 (2020). https://doi.org/10.1186/s13059-020-1930-x

Jingyun Luo, Chengcheng Wei, Haijun Liu, Shikun Cheng, Yingjie Xiao, Xiaqing Wang, Jianbing Yan, Jianxiao Liu, MaizeCUBIC: a comprehensive variation database for a maize synthetic population, Database, Volume 2020, 2020, baaa044, https://doi.org/10.1093/database/baaa044
}

\examples{
data("cubic_hmp")
}