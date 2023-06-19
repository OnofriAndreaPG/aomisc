\name{biplot.AMMIobj}
\alias{biplot.AMMIobject}
\alias{biplot.GGEobject}

\title{
Biplots for AMMI and GGE analyses of multi-environment genotype experiments
}
\description{
This method is reserved for the classes 'AMMIobject' and 'GGEobject' as obtained from the 'AMMI()' and 'GGE()' functions in the 'aomisc' package. It draws swift biplots of two types: type 1 plots PC1 against the average yield of genotypes across environments, while type 2 plots PC2 against PC1. The former type of biplot is reserved for 'AMMI' objects, while the second one is for AMMI and GGE objects.
}
\usage{
biplot.AMMIobj(obj, biplot = 1, xlim=NULL, ylim=NULL, elabels=NULL, 
               glabels=NULL, quad=FALSE, cexG=0.9, 
               cexE=0.9,
               xlab=NULL, ylab=NULL, font=1)
}
\arguments{
  \item{obj}{a 'AMMIobject' or 'GGEobject' objects}
  \item{biplot}{Numeric: either 1 or 2, to request on of the two available types of biplots (see description).}
  \item{xlim}{graphical parameter as in the 'plot' method}
  \item{ylim}{graphical parameter as in the 'plot' method}
  \item{elabels}{labels for the environments}
  \item{glabels}{labels for the genotypes}
  \item{quad}{logical. If TRUE, plots the axes}
  \item{cexG}{graphical parameter: the 'cex' for the genotype labels}
  \item{cexE}{graphical parameter: the 'cex' for the environment labels}
  \item{xlab}{graphical parameter as in the 'plot' method}
  \item{ylab}{graphical parameter as in the 'plot' method}
  \item{font}{graphical parameter as in the 'plot' method. Relating to the genotype labels}
}
\value{It only returns a graph}

\author{
Andrea Onofri
}
\examples{
data("WinterWheat")
tab <- with(WinterWheat, AMMI(Yield, Genotype, Year, Block, PC = 2))
biplot(tab)
}
