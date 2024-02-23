\name{anova.aovlist}
\alias{anova.aovlist}
\title{Prints an ANOVA table for an 'aovList' object}
\usage{
\method{anova}{aovlist}(object, ...)
}
\arguments{
  \item{object}{An object of class \code{"aovList"}}
  \item{...}{Other additional arguments}
}
\description{
  A wrapper for the 'summary' method for objects of class 'aovlist'.
}
\value{
  An object of class \code{"summary.aovlist"}. It is a list (one object 
  per error stratum) of ANOVA tables (class \code{"anova"}) with a row for
  each term in the model, plus one for \code{"Residuals"} if there
  are any.
}
\seealso{
  \code{\link{aov}}, \code{\link{summary}}, \code{\link{model.tables}},
  \code{\link{TukeyHSD}}
}
\examples{
# A split-plot design
rm(list=ls())
data(oats)
oats.aov <- aov(Y ~ N*V + Error(B/V), data = oats)
anova(oats.aov)
emmeans::emmeans(oats.aov, ~N)
}
\keyword{models}
\keyword{regression}