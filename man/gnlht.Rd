\name{gnlht}
\alias{gnlht}
\title{
Test nonlinear functions of model parameters
}
\description{
This function calculates linear/nonlinear functions of model parameters and returns an estimate with standard errors.
}
\usage{
gnlht(modelObj, funList, constList)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{modelObj}{a model object for which there are coef and vcov methods, and for which the named coefficient vector returned by coef is asymptotically normally distributed with asymptotic covariance matrix returned by vcov.}
  \item{funList}{a dataframe with one column, listing the  quoted strings that are the functions of the parameter estimates to be evaluated}
  \item{constList}{If necessary, a dataframe whose columns are the constants to be used in the calculations above. For each row of this dataframe, the functions above are evaluated}
  }
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
A data.frame
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
Andrea Onofri
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
# Load dataset
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
