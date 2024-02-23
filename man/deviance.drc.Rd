\name{deviance.drc}
\alias{deviance.drc}

\title{
Residual deviance for a non-linear regression fit.
}
\description{
Calculate the sum of squared residuals from a non-linear regression fit with the 'drm()' function in the 'drc()' package.
}
\usage{
\method{deviance}{drc}(object, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{a 'drc' object, for which the deviance is required}
  \item{...}{Other additional arguments, if necessary}
}

\value{
The value of the deviance extracted from the object object.
}
\author{
Andrea Onofri
}
\examples{
X <- c(1, 3, 5, 7, 9, 11, 13, 20)
Y <- c(8.22, 14.0, 17.2, 16.9, 19.2, 19.6, 19.4, 19.6)

# nls fit
model <- nls(Y ~ NLS.asymReg(X, init, m, plateau) )
deviance(model)

# drm fit
model2 <- drm(Y ~ X, fct = DRC.asymReg())
deviance(model2)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory (show via RShowDoc("KEYWORDS")):
% \keyword{ ~kwd1 }
% \keyword{ ~kwd2 }
% Use only one keyword per line.
% For non-standard keywords, use \concept instead of \keyword:
% \concept{ ~cpt1 }
% \concept{ ~cpt2 }
% Use only one concept per line.
