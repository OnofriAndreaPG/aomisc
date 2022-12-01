\name{expoDecay}
\alias{expoDecay.fun}
\alias{NLS.expoDecay}
\alias{NLS.expoDecay}

\title{Exponential decay function}

\description{
  These functions provide the exponential decay function (expoDecay), the exponential decay function with self-starter for the \code{\link{nls}} function and the exponential decay function with self-starter for the \code{\link{drm}} function in the drc package.  
}

\usage{
  expoDecay.fun(predictor, init, k)
  NLS.expoDecay(predictor, init, k)
  DRC.expoDecay(fixed = c(NA, NA), names = c("init", "k"), ...)
}

\arguments{
 \item{fixed}{numeric vector. Specifies which parameters are fixed and at what value they are fixed. NAs for parameter that are not fixed.}
 \item{names}{a vector of character strings giving the names of the parameters. The default is reasonable.}
}

\details{
  The exponential decay is given by the following function:
  \deqn{ f(x) = \textrm{init} \cdot \exp (- k x )\}}
  }
  
\value{
  expoDecay.fun and NLS.expoDecay return a numeric value, while DRC.expoDecay returns a list containing the nonlinear function, the self starter function and the parameter names.
}

\references{
Ratkowsky, DA (1990) Handbook of nonlinear regression models. New York (USA): Marcel Dekker Inc.
}

\author{Andrea Onofri}

\note{DRC.expoDecay is for use with the function \code{\link{drm}}.}

\seealso{}

\examples{
rm(list=ls())
devtools::load_all()
data("degradation")
model <- drm(Conc ~ Time, fct = DRC.expoDecay(),
             data = degradation)
summary(model)
model2 <- nls(Conc ~ NLS.expoDecay(Time, a, b), 
              data = degradation)
summary(model2)
}

\keyword{models}
\keyword{nonlinear}
\concept{}
