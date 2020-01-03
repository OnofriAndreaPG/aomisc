\name{gnlht}
\alias{gnlht}
\title{
Test nonlinear functions of model parameters
}
\description{
This function calculates linear/nonlinear functions of model parameters and returns an estimate with standard errors.
}
\usage{
gnlht(modelObj, funList, constList)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{modelObj}{a model object for which there are coef and vcov methods, and for which the named coefficient vector returned by coef is asymptotically normally distributed with asymptotic covariance matrix returned by vcov.}
  \item{funList}{a dataframe with one column, listing the  quoted strings that are the functions of the parameter estimates to be evaluated}
  \item{constList}{If necessary, a dataframe whose columns are the constants to be used in the calculations above. For each row of this dataframe, the functions above are evaluated}
  }
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
A data.frame
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
Andrea Onofri
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
# Load dataset
library(aomisc)
data(metamitron)

#Fit nls grouped model
modNlin <- nls(Conc ~ A[Tesi] * exp(-k[Tesi] * Tempo), 
               start=list(A=rep(100, 4), k=rep(0.06, 4)), 
               data=metamitron)
summary(modNlin)

#Make predictions
funList <- data.frame(form = c("A1 * exp(-k1 * Time)", "A2 * exp(-k2 * Time)"))
constList <- data.frame(Time = c(0, 7, 14, 21, 32, 42, 55, 67))
gnlht(modNlin, funList, constList)

#More complex predictions (2 constants)
funList <- data.frame(form = c("A1 * exp(-k1 * Time) * te", "A2 * exp(-k2 * Time) * te"))
constList <- data.frame(Time = c(0, 7, 14, 21, 32, 42, 55, 67), te = rep(0.5, 8))
gnlht(modNlin, funList, constList)

#More complex predictions (2 constants, but one function)
funList <- data.frame(form = c("A1 * exp(-k1 * Time) * te"))
constList <- data.frame(Time = c(0, 7, 14, 21, 32, 42, 55, 67), te = rep(0.5, 8))
gnlht(modNlin, funList, constList)

#Differences between parameters (no constants)
funList <- data.frame(form = c("k1 - k2", "k1 - k3", "k1 - k4"))
gnlht(modNlin, funList)

#Differences between parameters (no constants)
funList <- data.frame(form = c("k1 - k2"))
gnlht(modNlin, funList)

#Half-lives
funList <- data.frame(form = c("-log(0.5)/k1", "-log(0.5)/k2"))
gnlht(modNlin, funList)

#degradation times
funList <- data.frame(form = c("-log(tval)/k1"))
costList <- data.frame(tval = c(0.3, 0.5, 0.7))
gnlht(modNlin, funList, costList)

# Differences between half-lives
funList <- data.frame(form = c("-log(0.5)/k1 + log(0.5)/k2",
                               "-log(0.5)/k1 + log(0.5)/k3",
                               "-log(0.5)/k1 + log(0.5)/k4"))
gnlht(modNlin, funList)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
