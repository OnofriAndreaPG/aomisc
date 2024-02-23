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
\author{
Andrea Onofri
}

\examples{
fileName <- "https://www.casaonofri.it/_datasets/FGP_rape.csv"
library(aomisc)
dataset <- read.csv(fileName)
dataset[,1:5] <- lapply(dataset[,1:5], factor)
mod <- lm(FGP ~ Genotype * Run, data = dataset)
check <- check.hom(mod, Run)
check$aovtable
}
