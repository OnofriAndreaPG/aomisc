\name{linear}
\alias{linear.fun}
\alias{NLS.linear}
\alias{NLS.linearOrigin}
\alias{DRC.linear}

\title{Simple linear regression functions}

\description{
These functions provide the simple linear regression model (linear), the linear regression model with self-starter for the \code{\link{nls}} function (NLS.linear) and the simple linear regression function with self-starter for the \code{\link{drm}} function in the drc package (DRC.linear). For the 'nls' function, we also provide function and self starter for the simple linear regression through origin (NLS.linearOrigin). Obviously, fitting linear functions with nonlinear least square regression is sub-optimal, but it might be useful for comparing alternative models.
}

\usage{
  linear.fun(predictor, a, b)
  NLS.linear(predictor, a, b)
  NLS.linearOrigin(predictor, b)
  DRC.linear(fixed = c(NA, NA), names = c("a", "b"))
}

\arguments{
\item{predictor}{a numeric vector of values at which to evaluate the model.}
 \item{a}{numeric. The response when the predict is equal to 0.}
 \item{b}{numeric. The slope.}
 \item{fixed}{numeric vector. Specifies which parameters are fixed and at what value they are fixed. NAs for parameter that are not fixed.}
 \item{names}{a vector of character strings giving the names of the parameters. The default is reasonable.}
}

\details{
  The simple linear regression model is given by the following equation:
  \deqn{ f(x) = a +  b x }
  }
  
\value{
  linear.fun, NLS.linear and NLS.linearOrigin return a numeric value, while DRC.linear returns a list containing the nonlinear function, the self starter function and the parameter names.
}

\references{
Ratkowsky, DA (1990) Handbook of nonlinear regression models. New York (USA): Marcel Dekker Inc.

Onofri, A. (2020). A collection of self-starters for nonlinear regression in R. See: \url{https://www.statforbiology.com/2020/stat_nls_usefulfunctions/} 
}

\author{Andrea Onofri}

\examples{
# Simple linear regression
X <- seq(5, 50, 5)
Y <- 10 + 0.5*X + rnorm(10, 0, 0.5)

model1 <- nls(Y ~ NLS.linear(X, a, b))
model2 <- nls(Y ~ NLS.linearOrigin(X, b)) # force through origin
summary(model1); summary(model2)

model1 <- drm(Y ~ X, fct = DRC.linear())
model2 <- drm(Y ~ X, fct = DRC.linear(fixed = c(0, NA)))
summary(model1); summary(model2)
}

\keyword{models}
\keyword{nonlinear-regression}
