\name{expoGrowth}
\alias{expoGrowth.fun}
\alias{NLS.expoGrowth}
\alias{NLS.expoGrowth}

\title{Exponential growth function}

\description{
  These functions provide the exponential growth function (expoGrowth), the exponential growth function with self-starter for the \code{\link{nls}} function and the exponential growth function with self-starter for the \code{\link{drm}} function in the drc package.  
}

\usage{
  expoGrowth.fun(predictor, init, k)
  NLS.expoGrowth(predictor, init, k)
  DRC.expoGrowth(fixed = c(NA, NA), names = c("init", "k"), ...)
}

\arguments{
 \item{fixed}{numeric vector. Specifies which parameters are fixed and at what value they are fixed. NAs for parameter that are not fixed.}
 \item{names}{a vector of character strings giving the names of the parameters. The default is reasonable.}
}

\details{
  The exponential growth is given by the following function:
  \deqn{ f(x) = \textrm{init} \cdot \exp (k x )\}}
  }
  
\value{
  expoGrowth.fun and NLS.expoGrowth return a numeric value, while DRC.expoGrowth returns a list containing the nonlinear function, the self starter function and the parameter names.
}

\references{
Ratkowsky, DA (1990) Handbook of nonlinear regression models. New York (USA): Marcel Dekker Inc.
}

\author{Andrea Onofri}

\note{DRC.expoGrowth is for use with the function \code{\link{drm}}.}

\seealso{}

\examples{
rm(list=ls())
set.seed(1234)
X <- 1:20
Ye <- 3 * exp(0.1 * X)
Y <- Ye + rnorm(20, 0, 0.7)
dataset <- data.frame(X, Y)
rm(X, Y)

model <- drm(Y ~ X, fct = DRC.expoGrowth(),
             data = dataset)
model2 <- nls(Y ~ NLS.expoGrowth(X, a, b),
             data = dataset)
summary(model)
summary(model2)
}

\keyword{models}
\keyword{nonlinear}
\concept{}
