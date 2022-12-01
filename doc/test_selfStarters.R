library(devtools)
# install.packages("rlang")
install_github("OnofriAndreaPG/aomisc")
library(aomisc)

# Exponential growth ###################
devtools::load_all()
rm(list=ls())
set.seed(1234)
X <- 1:20
Ye <- 3 * exp(0.1 * X)
Y <- Ye + rnorm(20, 0, 0.7)
dataset <- data.frame(X, Y)
rm(X, Y)

model <- drm(Y ~ X, fct = DRC.expoGrowth(),
             data = dataset)
model2 <- nls(Y ~ NLS.expoGrowth(X, a, b),
             data = dataset)
summary(model)
summary(model2)
predict(model, se.fit = T)
predict(model2, se.fit = T) # Not working

# for nls, predictions can be done with
funList <- data.frame(form = c("a * exp(b * X)"))
constList <- data.frame(X = dataset$X)
gnlht(model2, funList, constList)

# with drc, we can also use the same method (slower)
gnlht(model, funList, constList, 
      parameterNames = c("a", "b"))

# Exponential decay ##########################
rm(list=ls())
devtools::load_all()
data("degradation")
model <- drm(Conc ~ Time, fct = DRC.expoDecay(),
             data = degradation)
summary(model)
model2 <- nls(Conc ~ NLS.expoDecay(Time, a, b), 
              data = degradation)
summary(model2)

predict(model, se.fit = T)
predict(model2, se.fit = T) # Not working

# for nls, predictions can be done with
funList <- data.frame(form = c("a * exp(-b * X)"))
constList <- data.frame(X = degradation$Time)
gnlht(model2, funList, constList)

# with drc, we can also use the same method (slower)
gnlht(model, funList, constList, 
      parameterNames = c("a", "b"))

# Half-life
devtools::load_all()
model <- drm(Conc ~ Time, fct = DRC.expoDecay(),
             data = degradation)
ED(model, respLev = c(10, 30, 50))

# Negative exponential function ##########
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



