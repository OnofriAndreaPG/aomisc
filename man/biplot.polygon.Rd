\name{biplot.polygon}
\alias{biplot.polygon}

\title{
Superimposes a which-won-where pattern to a GGE biplot
}
\description{
This function is reserved for the class 'GGEobject' as obtained from the 'GGE()' function in the 'aomisc' package. It needs a GGE biplot open in a graphic device and adds a polygon view, as defined in Yan and Tinker (2006).
}
\usage{
biplot.polygon(obj, vertex)
}
\arguments{
  \item{obj}{a 'GGEobject' objects}
  \item{vertex}{a vector indicating the number of row for the vertex genotypes, with reference to the table of genotypic scores in 'obj'}
}
\value{No value is returned}

\author{
Andrea Onofri
}
\examples{
data("WinterWheat")
tab <- with(WinterWheat, GGE(Yield, Genotype, Year, Block, PC = 2))
biplot(tab)
vertices <- c(1,2,7,3,5)
biplot.polygon(tab, vertex = vertices)
}
