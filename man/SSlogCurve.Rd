\name{logCurve}
\alias{logCurve.fun}
\alias{logCurveNI.fun}
\alias{NLS.logCurve}
\alias{NLS.logCurveNI}
\alias{DRC.logCurve}

\title{Logarithmic curve}

\description{
These functions provide the logarithmic model (logCurve) with self-starter for the \code{\link{nls}} function and for the \code{\link{drm}} function in the drc package.
}

\usage{
  logCurve.fun(predictor, a, b)
  NLS.logCurve(predictor, a, b)
  NLS.logCurveNI(predictor, b)
  DRC.logCurve(fixed = c(NA, NA), names = c("a", "b"))
}

\arguments{
\item{predictor}{a numeric vector of values at which to evaluate the model.}
\item{a}{model parameter}
 \item{b}{model parameter}
 \item{fixed}{numeric vector. Specifies which parameters are fixed and at what value they are fixed. NAs for parameter that are not fixed.}
 \item{names}{a vector of character strings giving the names of the parameters. The default is reasonable.}
}

\details{
  The logarithmic curve is given by the following function:
  \deqn{ f(x) = a + b \log (X) }
  This curve crosses the X axis at X = a. We can force it through the origin by setting a = 0; this is possible by setting 'fixed = c(=, NA), while, in the 'nls()' function, we need to use the NLS.logCurveNI()' function.
}
  
\value{
  logCurve.fun, NLS.logCurve and NLS.logCurveNI return a numeric value, while DRC.logCurve returns a list containing the nonlinear function, the self starter function and the parameter names.
}

\references{
Ratkowsky, DA (1990) Handbook of nonlinear regression models. New York (USA): Marcel Dekker Inc.

Onofri, A. (2020). A collection of self-starters for nonlinear regression in R. See: \url{https://www.statforbiology.com/2020/stat_nls_usefulfunctions/} 
}

\author{Andrea Onofri}

\note{DRC.logCurve() is for use with the function \code{\link{drm}}.}

\examples{
X <- c(1,2,4,5,7,12)
Y <- c(1.97, 2.32, 2.67, 2.71, 2.86, 3.09)

# lm fit
model <- lm(Y ~ log(X) )

# nls fit
model <- nls(Y ~ NLS.logCurve(X, a, b) )

# drm fit
model <- drm(Y ~ X, fct = DRC.logCurve() )
}

\keyword{models}
\keyword{nonlinear}
