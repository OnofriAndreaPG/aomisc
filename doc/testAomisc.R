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

