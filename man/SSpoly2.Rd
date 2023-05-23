\name{poly2}
\alias{poly2.fun}
\alias{NLS.poly2}
\alias{DRC.poly2}

\title{Simple polynomial regression functions}

\description{
These functions provide the simple polynomial (second order) regression model (poly2), the polynomial regression model with self-starter for the \code{\link{nls}} function (NLS.poly2) and the polynomial regression function with self-starter for the \code{\link{drm}} function in the drc package (DRC.poly2).  Fitting linear functions with nonlinear least square regression is sub-optimal, but it might be useful for comparing alternative models.
}

\usage{
  poly2.fun(predictor, a, b, c)
  NLS.poly2(predictor, a, b, c)
  DRC.linear(fixed = c(NA, NA, NA), names = c("a", "b", "c"), ...)
}

\arguments{
 \item{a}{numeric. The response when the predictor is equal to 0.}
 \item{b}{numeric. The slope at X = 0}
 \item{c}{numeric. Regression parameter}
 \item{fixed}{numeric vector. Specifies which parameters are fixed and at what value they are fixed. NAs for parameter that are not fixed.}
 \item{names}{a vector of character strings giving the names of the parameters. The default is reasonable.}
}

\details{
  The simple polynomial (second order) regression model is given by the following equation:
  \deqn{ f(x) = a +  b x + c x^2}
  }
  
\value{
  poly2.fun and NLS.poly2 return a numeric value, while DRC.poly2 returns a list containing the nonlinear function, the self starter function and the parameter names.
}

\references{
Ratkowsky, DA (1990) Handbook of nonlinear regression models. New York (USA): Marcel Dekker Inc.

Onofri, A. (2020). A collection of self-starters for nonlinear regression in R. See: \url{https://www.statforbiology.com/2020/stat_nls_usefulfunctions/} 
}

\author{Andrea Onofri}

\examples{
# Polynomial regression
X <- seq(5, 50, 5)
Y <- c(12.6, 74.1, 157.6, 225.5, 303.4, 462.8, 
       669.9, 805.3, 964.2, 1169)

model <- nls(Y ~ NLS.poly2(X, a, b, c))
summary(model)
model <- drm(Y ~ X, fct = DRC.poly2())
summary(model)
}

\keyword{models}
\keyword{nonlinear-regression}
