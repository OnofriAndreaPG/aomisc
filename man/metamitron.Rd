\docType{data}
\name{metamitron}
\alias{metamitron}
\title{Soil degradation for metamitro with different co-applied herbicides}
\format{A data frame with 96 observations on the following 3 variables.
\describe{
\item{Time}{a numeric vector, representing the time after the treatment}
 \item{Herbicide}{a factor with four levels, representing four herbicides} 
 \item{Conc}{a numeric vector, representing the concentrations} 
}}
\usage{data("metamitron")}
\source{
Vischetti, C., Marini, M., Businelli, M., Onofri, A., 1996. The effect of temperature and co-applied herbicides on the degradation rate of phenmedipham, chloridazon and metamitron in a clay loam soil in the laboratory, in: Re, A.D., Capri, E., Evans, S.P., Trevisan, M. (Eds.), “The Environmental Phate of Xenobiotics”, Proceedings X Symposium on Pesticide Chemistry, Piacenza. La Goliardica Pavese, Piacenza, pp. 287–294.
}
\description{
Herbicide degradation
}
\examples{
data(metamitron)
mod1 <- drm(Conc ~ Time, fct = DRC.expoDecay(),
           data = metamitron, curveid = Herbicide)
summary(mod1)
}

