\name{predictGS {G2P}}
\alias{predictGS}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Performing prediciton with fitted Models.
}
\description{
The function perform prediction from new inputs with fitted models, including RRBLUP, SVR, SVC, RFR, RFC, LASSO, SPLS and bigRR.
}
\usage{
predictGS(testMat, trainModel,predX = NULL,modelMethods = "SVC")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{testMat}{(numeric, matrix) a new set of genotypic data that need to predict. row represents sample and column represents SNP (feature).Genotypes generally be coded as {0,1,2};0 represents AA(homozygote),2 represents BB(homozygote) and 1 represents AB(heterozygote), the coding pattern must be consistent with training data. missing (NA) alleles are not allowed.}
  \item{trainModel}{(model) the fitted model. The object must be consistent with modelMethods.}
  \item{predX}{(numeric, matrix) a design matrix as fixed effects, the formation must be consistent with fitted model used.}
  \item{modelMethods}{(character) the name of fitted model, including RRBLUP, SVR, SVC, RFR, RFC, LASSO, SPLS and bigRR.}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
a list:\cr  
The prediction results of testing sets which predicted through already models.
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
## Load example data ##
data(cbuic)

## Fit RRBLUP model ##
RRBLUP_model <- GSReModel(markers = Markers[1:200,], pheVal = phenotype[1:200], modelMethods = "RRBLUP")

## Predict 1-20 subset of all example data with already RRBLUP model ## 
res <- predictGS(testMat = Markers[201:400,], trainModel = RRBLUP_model, modelMethods = "RRBLUP")
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
