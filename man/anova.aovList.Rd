\name{anova.aovList}
\alias{anova.aovList}
\title{Prints an ANOVA table for an 'aovList' object}
\usage{
\method{anova}{aovList}(object)
}
\arguments{
  \item{object}{An object of class \code{"aovList"}}
}
\description{
  A wrapper for the 'summary' method for objects of class 'aovList'.
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
# A spli-plot design
rm(list=ls())
data(oats)
oats.aov <- aov(Y ~ N*V + Error(B/V), data = oats)
anova(oats.aov)
TukeyHSD(oats.aov)
}
\keyword{models}
\keyword{regression}