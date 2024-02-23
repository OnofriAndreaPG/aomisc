\name{negExp}
\alias{negExp.fun}
\alias{negExpDist.fun}
\alias{NLS.negExp}
\alias{NLS.negExpDist}
\alias{DRC.negExp}
\alias{DRC.negExpDist}

\title{Negative exponential functions}

\description{
These functions provide the negative exponential model (negExp.fun) with the related self-starters for the \code{\link{nls}} function (NLS.negExp)  \code{\link{drm}} function in the 'drc' package (DRC.negExp) and the exponential cumulative distribution function (negExpDist.fun), with self-starters for both 'nls' (NLS.negExpDist) and 'drc' (DRC.negExpDist).
}

\usage{
  negExp.fun(predictor, a, c)
  negExpDist.fun(predictor, c)
  NLS.negExp(predictor, a, c)
  DRC.negExp(fixed = c(NA, NA), names = c("a", "c"))
  NLS.negExpDist(predictor, c)
  DRC.negExpDist(fixed = NA, names = c("c"))
}

\arguments{
\item{predictor}{a numeric vector of values at which to evaluate the model.}
\item{a}{a numeric parameter representing the higher asymptote}
\item{c}{a numeric parameter that is proportional to the relative rate of increase of the fitted function}
 \item{fixed}{numeric vector. Specifies which parameters are fixed and at what value they are fixed. NAs for parameter that are not fixed.}
 \item{names}{a vector of character strings giving the names of the parameters. The default is reasonable.}
}

\details{
  The negative exponential model is given by the following function:
  
  \deqn{ f(x) = \textrm{a} \{1 -  \exp [- \exp( c x )] \} }
  
  while the exponential CDF is obtained by setting a = 1:
  
  \deqn{ f(x) = 1 -  \exp [- \exp( c x )] }
  
  The ‘drc’ package contains also the function AR.2(), where c is replaced
  by e = 1/c. The ‘nlme’ package also contains an alternative
  parameterisation named 'SSasympOrig()', where c
  is replaced by ϕ3 = log(c).
  }
  
\value{
negExp.fun and negExpDist.fun return a numeric value, while the self-starters return a list containing the nonlinear function, the self starter function and the parameter names.
}

\references{
Ratkowsky, DA (1990) Handbook of nonlinear regression models. New York (USA): Marcel Dekker Inc.
Onofri, A. (2020). A collection of self-starters for nonlinear regression in R. See: \url{https://www.statforbiology.com/2020/stat_nls_usefulfunctions/} 
}

\author{Andrea Onofri}

\examples{
X <- c(1, 3, 5, 7, 9, 11, 13, 20)
Y <- c(4.5, 12.0, 16.1, 16.4, 18.9, 19.5, 19.3, 19.6)
model <- drm(Y ~ X, fct = DRC.negExp())
summary(model)
}

\keyword{models}
\keyword{nonlinear}