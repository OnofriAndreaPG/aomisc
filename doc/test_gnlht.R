# library(devtools)
# install_github("OnofriAndreaPG/aomisc")
rm(list = ls())
library(aomisc)
data(metamitron)

#Fit nls grouped model
modNlin <- nls(Conc ~ A[Herbicide] * exp(-k[Herbicide] * Time), 
               start=list(A=rep(100, 4), k=rep(0.06, 4)), 
               data=metamitron)
summary(modNlin)


# Compare parameters #######
funList <- list(~k1 - k2, ~k1 - k3, ~k1 - k4)
gnlht(modNlin, funList)

# Combine parameters
funList <- list(~ -log(0.5)/k1, ~ -log(0.5)/k2,
                ~ -log(0.5)/k3, ~ -log(0.5)/k4)
gnlht(modNlin, funList)

# Combine more flexibly
funList <- list(~ -log(prop)/k1, ~ -log(prop)/k2,
                ~ -log(prop)/k3, ~ -log(prop)/k4)
gnlht(modNlin, funList, const = data.frame(prop = 0.5))

funList <- list(~ -log(prop)/k1, ~ -log(prop)/k2,
                ~ -log(prop)/k3, ~ -log(prop)/k4)
gnlht(modNlin, funList, const = data.frame(prop = c(0.7, 0.5, 0.3)))

# Other possibilities
funList <- list(~ (k2 - k1)/(k1 * k2) * log(prop),
                ~ (k3 - k1)/(k1 * k3) * log(prop), 
                ~ (k4 - k1)/(k1 * k4) * log(prop))
gnlht(modNlin, funList, const = data.frame(prop = c(0.7, 0.5, 0.3)))

# Predictions
funList <- list(~ A1 * exp (- k1 * Time), ~ A2 * exp (- k2 * Time), 
                ~ A3 * exp (- k3 * Time), ~ A4 * exp (- k4 * Time))
propdF <- data.frame(Time = seq(0, 67, 1))
func <- funList
const <- propdF
pred <- gnlht(modNlin, funList, const = propdF)
head(pred)
tail(pred)

# Use with drm and nls
rm(list=ls())
library(aomisc)
data(speciesArea)
model <- drm(numSpecies ~ Area, fct = DRC.powerCurve(),
             data = speciesArea)
model2 <- nls(numSpecies ~ NLS.powerCurve(Area, a, b),
             data = speciesArea)

# First derivative
D(expression(a * X^b), "X")
## a * (X^(b - 1) * b)

# Second derivative
D(D(expression(a * X^b), "X"), "X")
## a * (X^((b - 1) - 1) * (b - 1) * b)

funList <- list(~ a * (X^(b - 1) * b), ~ a * (X^((b - 1) - 1) * (b - 1) * b))
propdF <- data.frame(X = seq(1, 5, 1))
pred <- gnlht(model, funList, const = propdF, parameterNames = c("a", "b"))
pred
pred2 <- gnlht(model2, funList, const = propdF, parameterNames = c("a", "b"))
pred2
# Compare parameters #######
funList <- list(~k1 - k2, ~k1 - k3, ~k1 - k4)
gnlht(modNlin, funList)

# Back transformation of means
library(emmeans)
data(mixture)
head(mixture)
mod <- lm(sqrt(Weight) ~ Treat, data = mixture)
med <- emmeans(mod, ~Treat)
coefs <- as.data.frame(med)[,2]
#parameterNames <- paste("b", 1:4, sep = "")
parameterNames <- as.data.frame(med)$Treat
pn <- parameterNames
# func <- as.list(paste("(1 -", pn, ")/", pn, sep = ""))
# func <- as.list(paste("exp(", pn, ") - 1", sep = ""))
# func <- as.list(paste("(1 -", pn, ")/", pn, sep = ""))
func <- as.list(paste(pn, "^2", sep = ""))
vcov. <- diag(as.data.frame(med)[,3]^2)
sigma <- vcov.
const <- NULL

medieBack <- gnlht(coefs, func,  vcov. = sigma,
                           parameterNames = pn, dfr = 12)
medieBack
emmeans(mod, ~Treat, transform = "response")

func <- as.list(paste(pn, "^2 + prova - pres", sep = ""))
vcov. <- diag(as.data.frame(med)[,3]^2)
sigma <- vcov.
const <- data.frame(prova = c(3, 4), pres = c(2, 3))
# const <- data.frame(prova = c(3, 4))
medieBack <- gnlht(coefs, func,  const = const, vcov. = sigma,
                           parameterNames = pn, dfr = 12)
medieBack

constants <- c(pres = 3, prova = 2)
car::deltaMethod(object=coefs, g = func[[1]], 
                      vcov.=vcMat, level = 0.95, 
                      constants = constants)
row.names(medieBack) <- as.character(as.data.frame(medie)[,1])
medieBack <- medieBack[,-1]

