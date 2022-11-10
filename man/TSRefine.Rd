\name{TSRefine {G2P}}
\alias{TSRefine}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Refinement of training set
}
\description{
Refining the training set for selective phenotyping to reduce cost.
}
\usage{
TSRefine(markers,candidates,test = NULL,ntosel,method = "PEVmean",
         npop = 100, nelite =5 ,mutprob = .8,niterations = 500,
         cores = 1,lambda = NULL,sequent = FALSE,visulization = F)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{markers}{(numeric, matrix) row represents sample and column represents SNP (feature).Genotypes generally be coded as {0,1,2};0 represents AA(homozygote),2 represents BB(homozygote) and 1 represents AB(heterozygote); missing (NA) alleles are not allowed.}
  \item{candidates}{(charcter array) names of condidates population.}
  \item{test}{(character, array) names of testing or predicted set, defult NULL, if NULL, execute reference-free refinement, else execute reference-based refinement.}
  \item{ntosel}{(integer) number of samples to select in refinement.}
  \item{method}{(character) alternativeal algorithms for training set refinement or selection, including "PEVmean", "CDmean", "sim", or "RD". Default is "PEVmean".}
  \item{npop}{(integer) genetic algorithm parameter, number of solutions at each iteration.}
  \item{nelite}{(integer) genetic algorithm parameter, number of solutions selected as elite parents which will generate the next set of solutions.}
  \item{mutprob}{(numeric) genetic algorithm parameter, probability of mutation for each generated soluation.}
  \item{niterations}{(integer) genetic algorithm parameter, number of iterations.}
  \item{cores}{(integer) number of CPU cores to use.}
  \item{lambda}{(numeric) scalar shrinkage parameter (λ>0).}
  \item{sequent}{(logical) phenotype selection parameter, if TRUE, perform sequential optimization, default FALSE.}
}
\details{

}
\value{
a list of selected samples with ID or index.
}
\references{
Sen, Ś., Satagopan, J. M., Broman, K. W., & Churchill, G. A. (2007). R/qtlDesign: inbred line cross experimental design. Mammalian genome, 18(2), 87-93.

Akdemir, D. (2017). STPGA: Selection of training populations with a genetic algorithm. BioRxiv, 111989.
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\link[STPGA]{GenAlgForSubsetSelection} \cr
\link[qtlDesign]{mma}
}
\examples{
## load example data ##
data(cubic)

## reference-free refinement ##
G <- Markers[1:300,]
TSR_SP <- TSRefine(markers = G, candidates = rownames(G),
              test = NULL, ntosel = 100,cores = 1,visulization = T)
              
## reference-based refinement (PEV-mean) ##
TSR_PEVmean <- TSRefine(markers = G, method = "PEVmean", candidates = rownames(G)[1:200],
                       test = rownames(Markers)[201:300], ntosel = 50,cores = 1,visulization = T)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory (show via RShowDoc("KEYWORDS")):
% \keyword{ ~kwd1 }
% \keyword{ ~kwd2 }
% Use only one keyword per line.
% For non-standard keywords, use \concept instead of \keyword:
% \concept{ ~cpt1 }
% \concept{ ~cpt2 }
% Use only one concept per line.
