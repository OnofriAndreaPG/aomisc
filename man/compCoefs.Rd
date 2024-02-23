\name{compCoefs}
\alias{compCoefs}
\title{
Pairwise comparisons of model parameters for nls objects
}
\description{
With models containing covariates and factors, we may be interested in fitting nonlinear regression models in a groupwise fashion, according to the levels of the experimental factor(s). This function permits pairwise comparison procedures of model parameters; parameters can be compared by means of either ratios or differences. It is based on the function 'compParm' in the package 'drc'.
}
\usage{
compCoefs(object, parameterNames, operator = "-", display = "pairwise",
          adjust = "holm", reversed = FALSE, level = 0.05, Letters = c(letters, LETTERS))
}
\arguments{
  \item{object}{ an object of class 'nls' }
  \item{parameterNames}{The name of parameters to compare. It must corresponds to the names of parameters in coefs(object)}
  \item{operator}{a character. If equal to "/" parameter ratios are compared. If equal to "-" (default) parameter differences are compared.}
  \item{display}{a character. If equal to "pairwise" (default), pairwise comparisons are displayed. Otherwise, a compact letter display is returned, according to the argument in 'level'.}
  \item{adjust}{A multiplicity adjustment method as in p.adjust. Defaults to "holm".}
  \item{reversed}{logical. Should the order of means/letters be decreasing? Defaults to FALSE, which means that the means and letters are in increasing order.}
  \item{level}{Protection level for compact letter display}
  \item{Letters}{Vector of distinct characters (or character strings) used to connect levels that are not significantly different. They should be recogizable when concatenated. The default behaviour is to use the small letters, followed by the capital letters. See help for 'multcompView::multcompLetters()' for futher detail)}
}
\value{returns a dataframe, containing the results}
\author{
Andrea Onofri
}

\examples{
data(metamitron)
modNlin <- nls(Conc ~ A[Herbicide] * exp(-k[Herbicide] * Time), 
               start=list(A=rep(100, 4), k=rep(0.06, 4)), 
               data=metamitron)
tab <- summary(modNlin)
pn <- c("k1", "k2", "k3", "k4")
compCoefs(modNlin, pn) # Holm multiplicity adjustment
compCoefs(modNlin, pn, adjust = "none")

# Displaying letters (P = 0.05 is default)
compCoefs(modNlin, pn, adjust = "none", 
        display = "cld")

# Calculate ratios
compCoefs(modNlin, pn, adjust = "none", 
        display = "pairwise", operator = "/")
}

