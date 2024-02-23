\name{CVA}
\alias{CVA}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Canonical variate analysis for multienvironment and multitrait genotype experiments}
\description{
This function performs canonical variate analysis as a descriptive visualisation tool. It is close to the 'lda()' function in the MASS package but it is not meant to be used for discriminant analyses.
}
\usage{
CVA(dataset, groups, scale = TRUE, constraint = 3)
}

\arguments{
  \item{dataset}{ dataset is a multidimensional matrix of observations }
  \item{groups}{ groups is a vector coding for groupings }
  \item{scale}{ whether the data needs to be standardised prior to analysis. Defaults to TRUE }
  \item{constraint}{ It is the type of scaling for eigenvectors, so that canonical variates have: 1 = unit within-group standard deviations (most common); 2 = unit total standard deviations; 3 = unit within group norms; 4 = unit total norms. It defaults to 3}
}
\details{
  More detail can be found in a blog page, at 'https://www.statforbiology.com/2023/stat_multivar_cva/'. Please, note that preliminary data transformations (e.g.: standardisation) are left to the user and must be performed prior to analyses (see example below).
}
\value{
  \item{TOT}{matrix of total variances-covariances}
  \item{B}{matrix of 'between-groups' variances-covariances}
  \item{W}{matrix of 'within-group' variances-covariances}
  \item{B/W}{matrix of W^{-1} B}
  \item{eigenvalues}{vector of eigenvalues}
  \item{eigenvectors}{matrix of eigenvectors}
  \item{proportion}{a vector containing the proportion of total discriminating ability captured by each canonical variate}
  \item{correlation}{vector of canonical correlations}
  \item{squared.canonical.correlation}{vector of squared canonical correlations}
  \item{coefficients}{matrix of canonical coefficients}
  \item{scores}{matrix of canonical scores}
  \item{centroids}{matrix of scores for centroids}
  \item{total.structure}{matrix of total canonical structure}
  \item{between.structure}{matrix of between-groups canonical structure}
  \item{within.structure}{matrix of within-groups canonical structure}
  \item{class.fun}{matrix of classifications functions}
  \item{class.val}{matrix of classification values}
  \item{within.structure}{matrix of within-groups canonical structure}
  \item{class}{vector of predicted classes}
}
\references{\href{https://www.statforbiology.com/2023/stat_multivar_cva/}{https://www.statforbiology.com/2023/stat_multivar_cva/}}
\author{Andrea Onofri}
\examples{
fileName <- "https://www.casaonofri.it/_datasets/WheatQuality4years.csv"
dataset <- read.csv(fileName)
dataset$Year <- factor(dataset$Year)
head(dataset)

# Standardise the data
groups <- dataset$Genotype
Z <- apply(dataset[,3:6], 2, scale, center = TRUE, scale = TRUE)
head(Z)

# Performs CVA
cvaobj <- CVA(Z, groups)
cvaobj
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{CVA}
\keyword{Multivariate statistics}
