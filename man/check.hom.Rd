\name{check.hom}
\alias{check.hom}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Check linear models for homoscedasticity
}
\description{
This function takes a linear model object as an argument and checks whether the residuals are homoscedastic, in relation to a stratification variable, that is given as an argument.
}
\usage{
check.hom(obj, var)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{obj}{a linear model object fitted with lm()}
  \item{var}{a vector containing a stratification variable}
}
\details{The function fits a gls with same structure as the input model, together with a heteroscedastic gls, where residuals are allowed to assume a different variance for each level of the stratification variable. The two models are compared by a LRT}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
\item{LRT}{the value of LRT}
\item{LRT}{the P-value of LRT}
\item{aovtable}{a summary table for the LRT}
\item{modHet}{the gls object containing the heteroscedastic fit}
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
Andrea Onofri
}

\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (x) 
{
  }
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
