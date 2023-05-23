\name{gnlht}
\alias{gnlht}
\title{
Linear/nonlinear functions of model parameters
}
\description{
This function calculates linear/nonlinear functions of model parameters and returns an estimate with standard errors.
}
\usage{
## S3 method for class 'numeric'
gnlht(coefs, func, const, vcov.,
            parameterNames, dfr)
## S3 method for class 'lm'
gnlht(model, func, const, ...)
## S3 method for class 'nls'
gnlht(model, func, const, ...)
## S3 method for class 'lme'
gnlht(model, func, const, ...)
## S3 method for class 'nlme'
gnlht(model, func, const, ...)
## S3 method for class 'drc'
gnlht(model, func, const, ...)
}
\arguments{
\item{coefs}{a numeric vector}
\item{vcov.}{a variance-covariance matrix, or a function to calculate it}
\item{parameterNames}{a character vector with namings for 'coefs'}
\item{dfr}{degrees of freedom}
\item{model}{a model object for which there are coef and vcov methods, and for which the named coefficient vector returned by coef is asymptotically normally distributed with asymptotic covariance matrix returned by vcov.}
  \item{model}{a model object for which there are coef and vcov methods, and for which the named coefficient vector returned by coef is asymptotically normally distributed with asymptotic covariance matrix returned by vcov.}
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
rm(list = ls())
library(aomisc)
data(metamitron)

#Fit nls grouped model
modNlin <- nls(Conc ~ A[Herbicide] * exp(-k[Herbicide] * Time), 
               start=list(A=rep(100, 4), k=rep(0.06, 4)), 
               data=metamitron)
summary(modNlin)


# Compare parameters #######
funList <- list(~k1 - k2, ~k1 - k3, ~k1 - k4)
gnlht(modNlin, funList)

# Combine parameters
funList <- list(~ -log(0.5)/k1, ~ -log(0.5)/k2,
                ~ -log(0.5)/k3, ~ -log(0.5)/k4)
gnlht(modNlin, funList)

# Combine more flexibly
funList <- list(~ -log(prop)/k1, ~ -log(prop)/k2,
                ~ -log(prop)/k3, ~ -log(prop)/k4)
gnlht(modNlin, funList, const = data.frame(prop = 0.5))

funList <- list(~ -log(prop)/k1, ~ -log(prop)/k2,
                ~ -log(prop)/k3, ~ -log(prop)/k4)
gnlht(modNlin, funList, const = data.frame(prop = c(0.7, 0.5, 0.3)))

# Other possibilities
funList <- list(~ (k2 - k1)/(k1 * k2) * log(prop),
                ~ (k3 - k1)/(k1 * k3) * log(prop), 
                ~ (k4 - k1)/(k1 * k4) * log(prop))
gnlht(modNlin, funList, const = data.frame(prop = c(0.7, 0.5, 0.3)))

# Predictions
funList <- list(~ A1 * exp (- k1 * Time), ~ A2 * exp (- k2 * Time), 
                ~ A3 * exp (- k3 * Time), ~ A4 * exp (- k4 * Time))
propdF <- data.frame(Time = seq(0, 67, 1))
func <- funList
const <- propdF
pred <- gnlht(modNlin, funList, const = propdF)
head(pred)
tail(pred)

# Back transformation of means
library(emmeans)
data(mixture)
head(mixture)
mod <- lm(sqrt(Weight) ~ Treat, data = mixture)
med <- emmeans(mod, ~Treat)
coefs <- as.data.frame(med)[,2]
#parameterNames <- paste("b", 1:4, sep = "")
parameterNames <- as.data.frame(med)$Treat
pn <- parameterNames
# func <- as.list(paste("(1 -", pn, ")/", pn, sep = ""))
# func <- as.list(paste("exp(", pn, ") - 1", sep = ""))
# func <- as.list(paste("(1 -", pn, ")/", pn, sep = ""))
func <- as.list(paste(pn, "^2", sep = ""))
vcov. <- diag(as.data.frame(med)[,3]^2)
sigma <- vcov.
const <- NULL

medieBack <- gnlht(coefs, func,  vcov. = sigma,
                           parameterNames = pn, dfr = 12)
medieBack
emmeans(mod, ~Treat, transform = "response")
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
