\name{GGE}
\alias{GGE}
\title{
GGE analysis for multienvironment genotype experiments
}
\description{
This function performs the GGE (Genotype plus Genotype by Environment interaction) analysis, according to Yan et al., 2000.}
\usage{
GGE(yield, genotype, environment, block, PC = 2)
}
\arguments{
  \item{yield}{a vector containing yield levels}
  \item{genotype}{a vector containing genotype codings}
  \item{environment}{a vector containing environment codings}
  \item{block}{a vector containing block codes for each environment}
  \item{PC}{the number of PCs that one wants to extract}
}
\value{ Returns a list of class 'GGEobject' with the following components
    \item{genotype_means}{The overall least squares genotype means}
    \item{environment_means}{The overall least squares environment means}
    \item{interaction_means}{The least squares means for the genotype by environment combinations}
    \item{ge_effect}{a two-way table of 'interaction effects'gge' effects}
    \item{additive_ANOVA}{an ANOVA table for the additive model}
    \item{GGE_summary}{a summary table for GGE analysis}
    \item{environment_scores}{a table of environment scores}
    \item{genotype_scores}{a table of genotype scores}
}
\references{
Yan, W., Hunt, L.A., Sheng, Q., Szlavnics, Z., 2000. Cultivar Evaluation and Mega-Environment Investigation Based on the GGE Biplot. Crop Science 40, 597–605.
}

\author{
Andrea Onofri
}

\examples{
data("WinterWheat")
tab <- with(WinterWheat, GGE(Yield, Genotype, Year, Block, PC = 2))
tab
}
\keyword{ GGE }
\keyword{ MET }

