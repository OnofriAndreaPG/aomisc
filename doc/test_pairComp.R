# library(devtools)
# install_github("OnofriAndreaPG/aomisc")
rm(list = ls())
library(aomisc)
data(metamitron)

#Fit nls grouped model
modNlin <- nls(Conc ~ A[Herbicide] * exp(-k[Herbicide] * Time), 
               start=list(A=rep(100, 4), k=rep(0.06, 4)), 
               data=metamitron)
tab <- summary(modNlin)
tab

# Compare parameters #######à
funList <- list(~k1 - k2, ~k1 - k3, ~k1 - k4)
gnlht(modNlin, funList)

# More efficiently
coefs <- tab$coef[5:8, 1]
SEs <- tab$coef[5:8, 2]
pn <- row.names(tab$coef)[5:8]
compVal(coefs, SEs, pn, df = tab$df[2]) # Holm multiplicity adjustment
compVal(coefs, SEs, pn, df = tab$df[2], adjust = "none")

# Displaying letters (P = 0.05 is default)
compVal(coefs, SEs, pn, df = tab$df[2], adjust = "none", 
        display = "cld")

# Calculate ratios
compVal(coefs, SEs, pn, df = tab$df[2], adjust = "none", 
        display = "pairwise", operator = "/")

# Another possibility (based on glht)
cp <- pairComp(coefs, SEs, pn, dfr = tab$df[2], adjust = "single-step")
cp$pairs

# Helper functions here (building pairwise comp mat)
tukeyMat(1:5, LETTERS[1:5])







