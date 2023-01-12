\name{asymReg}
\alias{asymReg.fun}
\alias{NLS.asymReg}
\alias{DRC.asymReg}
\alias{DRC.SSasymp}

\title{Asymptotic functions}

\description{
These functions provide the asymptotic regression model (asymReg), the asymptotic regression function with self-starter for the \code{\link{nls}} function and the asymptotic regression function with self-starter for the \code{\link{drm}} function in the drc package. The 'DRC.SSasymp()' function provides the self-starter for drc, to fit the same function as in the 'SSasymp()' function in the 'nlme' package. The asymptotic regression model is also known as the Mitscherlich law in agriculture and as the von Bertalanffy law in fisheries research.
}

\usage{
  asymReg.fun(predictor, init, m, plateau)
  NLS.asymReg(predictor, init, m, plateau)
  DRC.asymReg(fixed = c(NA, NA, NA), names = c("init", "m", "plateau"), ...)
  DRC.SSasymp(fixed = c(NA, NA, NA), names = c("Asym", "R0", "lrc"), ...)
}

\arguments{
 \item{fixed}{numeric vector. Specifies which parameters are fixed and at what value they are fixed. NAs for parameter that are not fixed.}
 \item{names}{a vector of character strings giving the names of the parameters. The default is reasonable.}
}

\details{
  The asymptotic model is given by the following function:
  \deqn{ f(x) = \textrm{plateau} - (\textrm{plateau} - \textrm{init}) \cdot \exp (- m x ) }
  A similar parameterisation where 'm' is replaced by exp(lrc) is provided in the 'nlme' package, with self-starter for the 'nls' function. Here, we provide the self-starter for the 'drm' function in the 'drc' package:
  \deqn{ f(x) = \textrm{Asym} + (\textrm{R0} - \textrm{Asym}) \cdot \exp [- exp(lrc) x ] }
  
  }
  
\value{
  asymReg.fun and NLS. asymReg return a numeric value, while DRC.asymReg and DRC.SSasymp return a list containing the nonlinear function, the self starter function and the parameter names.
}

\references{
Ratkowsky, DA (1990) Handbook of nonlinear regression models. New York (USA): Marcel Dekker Inc.

Onofri, A. (2020). A collection of self-starters for nonlinear regression in R. See: \url{https://www.statforbiology.com/2020/stat_nls_usefulfunctions/} 
}

\author{Andrea Onofri}

\note{DRC.asymReg and DRC.SSasymp are for use with the function \code{\link{drm}}.}

\seealso{}

\examples{
X <- c(1, 3, 5, 7, 9, 11, 13, 20)
Y <- c(8.22, 14.0, 17.2, 16.9, 19.2, 19.6, 19.4, 19.6)

# nls fit
model <- nls(Y ~ NLS.asymReg(X, init, m, plateau) )

# drm fit
model <- drm(Y ~ X, fct = DRC.asymReg())
summary(model)

model2 <- drm(Y ~ X, fct = DRC.SSasymp())
summary(model2)
}

\keyword{models}
\keyword{nonlinear}
\concept{}
