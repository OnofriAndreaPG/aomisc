#Negative exponential function
library(devtools)
install_github("OnofriAndreaPG/aomisc")
library(aomisc)
set.seed(1234)
X <- c(1, 3, 5, 7, 9, 11, 13, 20)
a <- 20; c <- 0.3

Ye <- negExp.fun(X, a, c)
epsilon <- rnorm(8, 0, 0.5)
Y <- Ye + epsilon

model <- drm(Y ~ X, fct = DRC.negExp())
summary(model)
plot(model, log="")

model <- nls(Y ~ NLS.negExp(X, a, c))
summary(model)
library(lattice)
plotnls(model, log="")

# Exponential decay
data("degradation")
model <- drm(Conc ~ Time, fct = DRC.expoDecay(),
             data = degradation)
summary(model)

model2 <- nls(Conc ~ NLS.expoDecay(Time, a, c), 
              data = degradation)
summary(model2)
              
