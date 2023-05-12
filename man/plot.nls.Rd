\name{plot.nls}

\alias{plot.nls}

\title{Plotting diagnostics for an \code{nls} object}

\description{
  Three plots (selectable by the argument 'which') are currently available: a plot of residuals against fitted values (which = 1), a Normal Q-Q plot (which = 2) and a plot of predicted against expected (line) and observed (symbols). By default, the first graph is provided. As for the third graph, we can either plot all the data (type= "all") or the group means (type = "means"; the default)
}

\usage{
  
  \method{plot}{nls}(fm, type = "means", 
                     xlim = NULL, res = 100, 
                     which = 3, ...)
}

\arguments{
  \item{fm}{an object of class 'nls'}
  \item{type}{it can be either "means" or "all". In the first case, the group means are plotted for the third graph. It is iìonly used when which = 3}
  \item{xlim}{The limits for the x-axis (x1, x2)}
  \item{res}{For 'which = 3', it is the resolution of the fitted line}
  \item{which}{The type of graph: can be 1, 2 or 3 (see description)}
  \item{...}{additional graphical arguments}
}

\details{
  It mimicks the behaviour of the function plot.lm()
}

\author{
  Andrea Onofri 
}

\examples{
data(degradation)
mod <- nls(Conc ~ NLS.expoDecay(Time, a, k),
               data = degradation)
plot(mod, which = 3)
}

