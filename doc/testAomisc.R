

# gnlht #######################################################
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


# compVal #######################################
rm(list = ls())
data("metamitron")
head(metamitron)
library(dplyr)
dataset <- metamitron %>% 
  mutate(Time = factor(Time))

mod <- lm(Conc ~ Time*Herbicide, data = dataset)

library(emmeans)
medie <- emmeans(mod, ~Time:Herbicide)
contrast(medie, method = "pairwise", adjust = "holm")
multcomp::cld(medie, adjust = "holm", Letters = letters)
contrast(medie, method = "pairwise", adjust = "none")
multcomp::cld(medie, adjust = "none", Letters = letters)

obj <- as.data.frame(medie)[,3]
SE <- as.data.frame(medie)[,4]
parNam <- paste(as.data.frame(medie)[,1], as.data.frame(medie)[,2], sep = ":")
tab <- compVal(obj, parNam, SE, adjust = "none", method = "pairwise",
               df.residual = 64)
tab

tab <- compVal(obj, parNam, SE, adjust = "none", method = "cld",
               df.residual = 64, decreasing = F)
tab



