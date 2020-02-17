library(devtools)
install.packages("rlang")
install_github("OnofriAndreaPG/aomisc")
library(aomisc)


#Negative exponential function
rm(list=ls())
set.seed(1234)
X <- c(1, 3, 5, 7, 9, 11, 13, 20)
a <- 20; c <- 0.3
Ye <- negExp.fun(X, a, c)
epsilon <- rnorm(8, 0, 0.5)
Y <- Ye + epsilon

model <- drm(Y ~ X, fct = DRC.negExp())
summary(model)
plot(model, log="")
boxcox(model)

mod <- lm(Y ~ X)
p <- boxcox(mod)
p

model <- nls(Y ~ NLS.negExp(X, a, c))
summary(model)
plotnls(model, log="")
boxcox(model)

dataset <- data.frame(X, Y)
rm(X, Y, Ye, epsilon, a, c)
model <- nls(Y ~ NLS.negExp(X, a, c), data = dataset)
boxcox(model)

# Exponential decay
rm(list=ls())
data("degradation")
model <- drm(Conc ~ Time, fct = DRC.expoDecay(),
             data = degradation)
summary(model)
boxcox(model)
model2 <- nls(Conc ~ NLS.expoDecay(Time, a, c), 
              data = degradation)
summary(model2)
boxcox(model2)              

model2 <- nls(log(Conc) ~ log(NLS.expoDecay(Time, a, c)), 
              data = degradation, start = coef(model2))

model2 <- nls(Conc^0.5 ~ NLS.expoDecay(Time, a, c)^0.5, 
              data = degradation, start = coef(model2))

summary(model2)
plotnls(model2)
bcSummary(p)
str(p)
p$m$Rmat()
str <- model2$m$formula()[[2]]
str <- paste( "I(", str, "^0.5 - 1)/0.5", sep ="")
model2$m$formula()[[2]] <- formula(str)
newdata <- eval(model2$data)
p <- nls(formula(str ~ ((NLS.expoDecay(Time, a, c))^0.5 - 1)/0.5), 
            start = coef(model2), trace = FALSE, data = newdata)
summary(p)
library(lattice)
plotnls(p)
# gnlht
library(aomisc)
data(metamitron)

#Fit nls grouped model
modNlin <- nls(Conc ~ A[Tesi] * exp(-k[Tesi] * Tempo), 
               start=list(A=rep(100, 4), k=rep(0.06, 4)), 
               data=metamitron)
summary(modNlin)

#Make predictions
funList <- data.frame(form = c("A1 * exp(-k1 * Time)", "A2 * exp(-k2 * Time)"))
constList <- data.frame(Time = c(0, 7, 14, 21, 32, 42, 55, 67))
gnlht(modNlin, funList, constList)

#More complex predictions (2 constants)
funList <- data.frame(form = c("A1 * exp(-k1 * Time) * te", "A2 * exp(-k2 * Time) * te"))
constList <- data.frame(Time = c(0, 7, 14, 21, 32, 42, 55, 67), te = rep(0.5, 8))
gnlht(modNlin, funList, constList)

#More complex predictions (2 constants, but one function)
funList <- data.frame(form = c("A1 * exp(-k1 * Time) * te"))
constList <- data.frame(Time = c(0, 7, 14, 21, 32, 42, 55, 67), te = rep(0.5, 8))
gnlht(modNlin, funList, constList)

#Differences between parameters (no constants)
funList <- data.frame(form = c("k1 - k2", "k1 - k3", "k1 - k4"))
gnlht(modNlin, funList)

#Differences between parameters (no constants)
funList <- data.frame(form = c("k1 - k2"))
gnlht(modNlin, funList)

#Half-lives
funList <- data.frame(form = c("-log(0.5)/k1", "-log(0.5)/k2"))
gnlht(modNlin, funList)

#degradation times
funList <- data.frame(form = c("-log(tval)/k1"))
costList <- data.frame(tval = c(0.3, 0.5, 0.7))
gnlht(modNlin, funList, costList)

# Differences between half-lives
funList <- data.frame(form = c("-log(0.5)/k1 + log(0.5)/k2",
                               "-log(0.5)/k1 + log(0.5)/k3",
                               "-log(0.5)/k1 + log(0.5)/k4"))
gnlht(modNlin, funList)

