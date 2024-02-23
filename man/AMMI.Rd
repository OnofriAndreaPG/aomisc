\name{AMMI}
\alias{AMMI}
\title{
AMMI analysis for multienvironment genotype experiments
}
\description{
This function performs the AMMI (Addittive Main effects Multiplicative Interaction) analysis, according to Zobel et al (1988). The function has been described in Onofri and Ciriciofolo (2007).
}
\usage{
AMMI(yield, genotype, environment, block = NULL, PC = 2,
     MSE = NULL, dfr = NULL)
}
\arguments{
  \item{yield}{a vector containing yield levels}
  \item{genotype}{a vector containing genotype codings}
  \item{environment}{a vector containing environment codings}
  \item{block}{a vector containing block codes for each environment}
  \item{PC}{the number of PCs that one wants to extract}
  \item{MSE}{Mean Squared Error}
  \item{dfr}{Residual Degrees of Freedom}
}
\value{ Returns a list of class 'AMMIobject' with the following components
    \item{genotype_means}{The overall least squares genotype means}
    \item{environment_means}{The overall least squares environment means}
    \item{interaction_means}{The least squares means for the genotype by environment combinations}
    \item{interaction_effect}{a two-way table of interaction effects}
    \item{additive_ANOVA}{an ANOVA table for the additive model}
    \item{mult_Interaction}{an ANOVA table for multiplicative model}
    \item{MSE}{Mean Square Error}
    \item{dfr}{Degrees of freedom for the MSE}
    \item{environment_scores}{a table of environment scores}
    \item{genotype_scores}{a table of genotype scores}
    \item{stability}{AMMI stability value (ASV; Mohammadi and Amri, 2008)}
}
\references{
Mohammadi, R., Amri, A., 2008. Comparison of parametric and non-parametric methods for selecting stable and adapted durum wheat genotypes in varibale environments. Euphytica 159, 419–432.

Onofri, A., Ciriciofolo, E., 2007. Using R to perform the AMMI analysis on agriculture variety trials. R NEWS 7, 14–19.

Zobel, R. W., Wright, M.J., and Gauch, H. G., 1988. Statistical analysis of a yield trial. Agronomy Journal, 388-393.
}

\author{
Andrea Onofri
}

\examples{
data("WinterWheat")
tab <- with(WinterWheat, AMMI(Yield, Genotype, Year, Block, PC = 2))
tab
}
\keyword{ AMMI }
\keyword{ MET }

