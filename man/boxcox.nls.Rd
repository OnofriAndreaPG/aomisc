\name{boxcox.nls}
\alias{boxcox.nls}
\alias{bcSummary}

\title{Transform-both-sides Box-Cox transformation}

\description{
  Finds the optimal Box-Cox transformation for non-linear regression models.
}

\usage{
  boxcox.nls(object, lambda = seq(-2, 2, 1/10), plotit = TRUE, start, eps = 1/50, bcAdd = 0, level = 0.95, 
  xlab = expression(lambda), ylab = "log-likelihood", ...)

  bcSummary(object)
}

\arguments{
  \item{object}{object of class \code{\link{nls}}. For \code{bcSummary} the \code{nls} fit should have been
  obtained using \code{boxcox.nls}}
  \item{lambda}{numeric vector of lambda values; the default is (-2, 2) in steps of 0.1.}
  \item{plotit}{logical which controls whether the result should be plotted.}
  \item{start}{a list of starting values (optional).}
  \item{eps}{numeric value: the tolerance for lambda = 0; defaults to 0.02.}
  \item{bcAdd}{numeric value specifying the constant to be added on both sides prior to Box-Cox transformation. 
  The default is 0.}
  \item{level}{numeric value: the confidence level required.}
  \item{xlab}{character string: the label on the x axis, defaults to "lambda".}
  \item{ylab}{character string: the label on the y axis, defaults to "log-likelihood".}
  \item{\dots}{additional graphical parameters.}
}

\details{
  \code{boxcox.nls} is very similar to the \code{\link[MASS]{boxcox}} in its arguments.

  The optimal lambda value is determined using a profile likelihood approach: 
  For each lambda value the non-linear regression model is fitted and the lambda value resulting in thre largest
  value of the log likelihood function is picked.
  
  If a self starter model was used in the model fit, then gradient information will be used in the
  profiling.
}

\value{
  An object of class \code{\link{nls}} (returned invisibly). 
  If plotit = TRUE a plot of loglik vs lambda is shown indicating a confidence interval (by default 95\%) about 
  the optimal lambda value.
}

\references{
  Carroll, R. J. and Ruppert, D. (1988) \emph{Transformation and Weighting in Regression}, 
  New York: Chapman and Hall (Chapter 4).
}

\author{Christian Ritz, modified by Andrea Onofri}

%\note{}

\seealso{For linear regression the analogue is \code{\link[MASS]{boxcox}}.}

\examples{

## Fitting log-logistic model without transformation
ryegrass.m1<-nls(rootl~c+(d-c)/(1+(conc/e)^b),data=ryegrass, start=list(b=1,c=0,d=8,e=3))
summary(ryegrass.m1)

## Fitting the same model with optimal Box-Cox transformation
ryegrass.m2 <- boxcox(ryegrass.m1)
summary(ryegrass.m2)
bcSummary(ryegrass.m2)

## Fitting the Michaelis-Menten model without self starter
L.minor.m1 <- nls(rate ~ Vm*conc/(K+conc), data = L.minor, start = list(K=20, Vm=120))
L.minor.m2 <- boxcox(L.minor.m1)
bcSummary(L.minor.m2)

## Fitting the Michaelis-Menten model with self starter
L.minor.m3 <- nls(rate ~ SSmicmen(conc, Vm, K), data = L.minor)
L.minor.m4 <- boxcox(L.minor.m3)
bcSummary(L.minor.m4)

}

\keyword{models}
\keyword{nonlinear}

\concept{Box-Cox transform-both-sides transformation}
