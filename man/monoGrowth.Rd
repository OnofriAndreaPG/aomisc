\name{monoGrowth}

\alias{monoGrowth}

\title{Monomolecular growth model}

\description{
  Monomolecular growth model, with biologically meaningful parameters
}

\usage{
monoGrowth(fixed = c(NA, NA, NA), names = c("init", "m", "plateau"))
}

\arguments{
  \item{fixed}{numeric vector. Specifies which parameters are fixed and at what value they are fixed. 
  NAs for parameter that are not fixed.}
  \item{names}{vector of character strings giving the names of the parameters (should not contain ":").
  The default parameter names are: init, m, plateau.}
}

\details{
  The monomulecular growth model is an asymptotically increasing curve, that has been reparameterised to include some biologically 
  meaningful parameters, such as initial and final plant weight and initial relative growth rate. The mean function is:
  
\deqn{ f(x) = plateau - (plateau - init) * exp (- m * x) }
  
  The parameter init is the initial plant weight (attained at x=0), the parameter plateau is the final plant weight 
  reached for x going to infinity and the parameter m is the initial value of relative growth rate (when x=0).
  Thus the curve is monotonously increasing in x.
}

\value{
  A list of class \code{drcMean}, containing the mean function, the self starter function,
  the parameter names.
}

\references{Roderick Hunt, 1982. Plant Growth Curves. Edward Arnold Publisher, Great Britain, 248 pp}

\author{Andrea Onofri}

\note{Growth functions are generally fitted on log-transformed weight data, which equals to setting bc parameter to 0}

\seealso{}

\examples{

## Fitting a monomolacolar growth curve
load("tomato.Rdata")
model.3 <- drm(weight ~ DAT, data=tom2003, fct=monoGrowth())
plot(model.3, log="")
summary(model.3)
}
\keyword{growth}
\keyword{models}
\keyword{nonlinear}
