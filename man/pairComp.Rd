\name{pairComp}
\alias{pairComp}
\title{
Pairwise comparisons between the elements of a vector
}
\description{
This function provides pairwise comparisons between the element of a vector, as long as a variance-covariance matrix is also provided
}
\usage{
pairComp(parm, vcovMat, nams = NULL, dfr = NULL, adjust = "none",
         level = 0.05, Letters = c(letters, LETTERS, "."), 
         reversed = FALSE)
}
\arguments{
  \item{parm}{
A (possibly named) vector of estimates
}
  \item{vcovMat}{
Variance-covariance matrix for the estimates
}
  \item{nams}{
A character vector of parameter names (optional). If it is not provided and if 'parm' is not a named vectors, numbers '1:length(parm)' are used.
}
  \item{dfr}{
An optional number of residual degrees of freedom (defaults to Inf)
}
  \item{adjust}{
A multiplicity adjustment method as in 'multcomp'. Defaults to "none".
}
  \item{level}{
Protection level for compact letter display
}
  \item{Letters}{
Vector of distinct characters (or character strings) used to connect levels that are not significantly different. They should be recogizable when concatenated. The default behaviour is to use the small letters, followed by the capital letters. See help for 'multcompView::multcompLetters()' for futher detail)
}
  \item{reversed}{logical. Should the order of means/letters be decreasing? Defaults to FALSE, which means that the means and letters are in increasing order)
}
}
\value{
Returns a list with the following elements
\item{pairs}{A dataframe of pairwise comparisons}
\item{Letters}{A dataframe with compact letter display}
}
\references{
Onofri A. (2020) The broken bridge between biologists and statisticians: a blog and R package, Statforbiology, IT, web: https://www.statforbiology.com
}
\author{
Andrea Onofri
}

\examples{
# library(devtools)
# install_github("OnofriAndreaPG/aomisc")
library(aomisc)
data(metamitron)

#Fit nls grouped model
modNlin <- nls(Conc ~ A[Herbicide] * exp(-k[Herbicide] * Time), 
               start=list(A=rep(100, 4), k=rep(0.06, 4)), 
               data=metamitron)
summary(modNlin)

# Retreive infos and make comparisons
coefs <- coef(modNlin)[5:8] 
vcovMat <- vcov(modNlin)[5:8, 5:8]
cp <- pairComp(coefs, vcovMat, dfr = tab$df[2], 
               adjust = "none", reversed = F)
cp$Letters
cp$pairs
}
