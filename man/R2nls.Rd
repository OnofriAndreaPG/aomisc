\name{R2nls}
\alias{R2nls}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Goodness of fit for nonlinear regression
}
\description{
This function calculates measures  of goodness of fit for nonlinear regression. It works with both 'nls' and 'drc' objects
}
\usage{
R2nls(object)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object }{A nonlinear regression fit object. It can be either a 'nls' fit or 'drm' fit.
  }
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
  A list with the following slots:
  \item{R2 }{Traditional coefficient of determination, calculated as the ratio of model SS to total SS. Formula as in Schabenberger and Pierce, 5.23, pag 211.}
  \item{PseudoR2 }{Pseudo-R2, more useful for nonlinear regression with no-intercept-models. Formula Formula as in Schabenberger and Pierce, 5.24, pag 212.}
 \item{R2adj }{Adjusted R2, similar to R2 above, but penalised for higher number of parameters.}
 \item{MSE }{Mean Squared Error}
 \item{RMSE }{ Root Means Squared Error}
 \item{RRMSE }{Relative Root Means Squared Error}
 
}
\references{
Schabenberger, O., Pierce, F.J., 2002. Contemporary statistical models for the plant and soil sciences. Taylor & Francis, CRC Press, Books.

}
\author{
Andrea Onofri
}
\examples{
data(beetGrowth)
mod3 <- nls(weightInf ~ NLS.L3(DAE, b, c, d), data = beetGrowth)
R2nls(mod3)

mod4 <- drm(weightInf ~ DAE, fct = L.3(), data = beetGrowth)
R2nls(mod4)
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
