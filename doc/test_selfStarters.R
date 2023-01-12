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
predict(model, se.fit = T, interval = "confidence")
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

# Asymptotic regression #####################
X <- c(1, 3, 5, 7, 9, 11, 13, 20)
Y <- c(8.22, 14.0, 17.2, 16.9, 19.2, 19.6, 19.4, 19.6)

# nls fit
model1 <- nls(Y ~ NLS.asymReg(X, init, m, plateau) )
# for nls, predictions can be done with
funList <- data.frame(form = c("c - (c - a) * exp(- b * X)"))
constList <- data.frame(X = X)
gnlht(model1, funList, constList,
      parameterNames = c("a", "b", "c"))

# drm fit
model2 <- drm(Y ~ X, fct = DRC.asymReg())
plot(model2, log="", main = "Asymptotic regression", 
     ylim = c(0,25), xlim = c(0,20))
summary(model2)
predict(model2, se.fit = T)
# with drc, we can also use the same method (slower)
gnlht(model2, funList, constList, 
      parameterNames = c("a", "b", "c"))

model3 <- drm(Y ~ X, fct = DRC.SSasymp())
plot(model3, log="", main = "Asymptotic regression - 2", 
     ylim = c(0,25), xlim = c(0,20))
summary(model3)
predict(model3, se.fit = T)
funList <- data.frame(form = c("a + (b - a) * exp(- exp(c) * X)"))
gnlht(model3, funList, constList, 
      parameterNames = c("a", "b", "c"))

# Negative exponential function ##########
rm(list=ls())
devtools::load_all()
set.seed(1234)
X <- c(1, 3, 5, 7, 9, 11, 13, 20)
a <- 20; c <- 0.3
Ye <- negExp.fun(X, a, c)
epsilon <- rnorm(8, 0, 0.5)
Y <- Ye + epsilon

model <- drm(Y ~ X, fct = DRC.negExp())
summary(model)
plot(model, log="")
predict(model, se.fit = T)
funList <- data.frame(form = c("a * (1 - exp(-b * X))"))
constList <- data.frame(X = X)
gnlht(model, funList, constList, 
      parameterNames = c("a", "b"))

# Alternative
model2 <- drm(Y ~ X, fct = DRC.asymReg(fixed = c(0, NA, NA)))
summary(model2)
predict(model2, se.fit = T)

# Negative exponential distribution function ##########
rm(list=ls())
devtools::load_all()
set.seed(1234)
X <- c(1, 3, 5, 7, 9, 11, 13, 20)
a <- 1; c <- 0.3
Ye <- negExp.fun(X, a, c)
epsilon <- rnorm(8, 0, 0.03)
Y <- Ye + epsilon

model <- drm(Y ~ X, fct = DRC.negExp(fixed = c(NA, NA)))
summary(model)
plot(model, log="")
predict(model, se.fit = T)

model$fct$deriv1(X, matrix(coef(model), ncol=1))
model$fct$derivx(X, matrix(coef(model), ncol=1))

funList <- data.frame(form = c("a * (1 - exp(-b * X))"))
constList <- data.frame(X = X)
gnlht(model, funList, constList, 
      parameterNames = c("a", "b"))

# Alternative
model2 <- drm(Y ~ X, fct = DRC.asymReg(fixed = c(0, NA, NA)))
summary(model2)
predict(model2, se.fit = T)

# Logarithmic curve ###########################
rm(list = ls())
devtools::load_all()
X <- c(1,2,4,5,7,12)
Y <- c(1.97, 2.32, 2.67, 2.71, 2.86, 3.09)

# lm fit
model <- lm(Y ~ log(X) )

# nls fit
model <- nls(Y ~ NLS.logCurve(X, a, b) )
summary(model)

# drm fit
model <- drm(Y ~ X, fct = DRC.logCurve() )
summary(model)
predict(model, se.fit = T)
funList <- data.frame(form = c("a + b * log(X)"))
constList <- data.frame(X = X)
gnlht(model, funList, constList, 
      parameterNames = c("a", "b"))

# Alternative
model2 <- drm(Y ~ X, fct = DRC.asymReg(fixed = c(0, NA, NA)))
summary(model2)
predict(model2, se.fit = T)

# force through the origin (bad fit!)
dataset <- data.frame(X, Y)
model <- nls(Y ~ NLS.logCurveNI(X, b) )
summary(model)

# drm fit
model <- drm(Y ~ X, fct = DRC.logCurve(fixed = c(0, NA)) )
summary(model)
plot(model, log = "")


