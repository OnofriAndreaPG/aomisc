fileName <- "https://www.casaonofri.it/_datasets/degradation.csv"
degradation <- read.csv(fileName)
head(degradation, 10)
modNlin <- nls(Conc ~ A*exp(-k*Time), 
               start=list(A=100, k=0.05), 
               data=degradation)
summary(modNlin)
plot(modNlin, which = 1)
