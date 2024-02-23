\name{BC.3.fun}
\alias{BC.3.fun}
\alias{BC.4.fun}
\alias{BC.5.fun}

\title{The Brain-Cousens hormesis models}

\description{
These functions provide the Brain-Cousens models for describing u-shaped hormesis. They can be used with 'nls()' or 'nlme()' or 'gnls()' and they do not contain self-starting routines.
}

\usage{
  BC.5.fun(predictor, b, c, d, e, f)
  BC.4.fun(predictor, b, d, e, f)
  BC.3.fun(predictor, b, e, f)
}

\arguments{
  \item{predictor}{a numeric vector of values at which to evaluate the model}
\item{b}{model parameter: slope near to the inflection point}
\item{c}{model parameter: lower asymptote}
\item{d}{model parameter: higher asymptote}
\item{e}{model parameter: dose at inflection point}
\item{f}{model parameter for hormesis}
}

\details{
  The model functions for the Brain-Cousens model (Brain and Cousens, 1989) is
  \deqn{ f(x, b,c,d,e,f) = c + \frac{d-c+fx}{1+\exp(b(\log(x)-\log(e)))}},
  \deqn{ f(x, b,d,e,f) = c + \frac{d+fx}{1+\exp(b(\log(x)-\log(e)))}},
  \deqn{ f(x, b,e,f) = c + \frac{1+fx}{1+\exp(b(\log(x)-\log(e)))}},
These functions do not provide a selfStarting routine, but you can obtain starting values from a drc fit, by using the corresponding functions BC.4() and BC.5().
}
  
\value{
These functions return a numeric value
}

\references{
Ratkowsky, DA (1990) Handbook of nonlinear regression models. New York (USA): Marcel Dekker Inc.
Onofri, A. (2020). A collection of self-starters for nonlinear regression in R. See: \url{https://www.statforbiology.com/2020/stat_nls_usefulfunctions/} 
Brain, P. and Cousens, R. (1989) An equation to describe dose responses 
  where there is stimulation of growth at low doses, 
  \emph{Weed Research}, \bold{29}, 93--96.
}

\author{Andrea Onofri}


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
