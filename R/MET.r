summary.met <- function(Yield, Var, Env){
# Questa funzione calcola le statistiche descrittive per prove varietali
# ripetute in molti ambienti, anche sbilanciate. ATTENZIONE: lavora sulle 
# medie ambientali delle varieta (non sulle repliche).
# Data: 15/04/2009

#Statistiche descrittive
result <- data.frame("Count"=tapply(Yield, Var, length))
result <- cbind(result, "Media"=tapply(Yield, Var, mean))
result <- cbind(result, "Max"=tapply(Yield, Var, max))
result <- cbind(result, "Min"=tapply(Yield, Var, min))
result <- cbind(result, "SD"=tapply(Yield, Var, sd))

#lsMeans
options(contrasts=c('contr.sum','contr.poly'))
v <- factor(Var); e <- factor(Env)
model <- aov(Yield ~ v + e)
model.tables(model, type="means")
NumVar <- length(levels(v))
NumEnv <- length(levels(e))
Int <- coef(model)[1]
EffVar <- c(coef(model)[2:NumVar], -sum(coef(model)[2:NumVar])) 
EffEnv <- c(coef(model)[(1+NumVar):(1+NumVar-1+NumEnv-1)], -sum(coef(model)[(1+NumVar):(1+NumVar-1+NumEnv-1)]))
MedieVar <- Int + EffVar
result <- cbind(result, "LSmeans"=MedieVar)
result <- result[order(-result$Count), ]


#LSMeans Mixed
library(nlme)
options(contrasts=c('contr.sum','contr.poly'))
model.mix <- lme(Yield ~ v, random = ~ 1|e)
BLUEs <- fixef(model.mix)[1] + c(fixef(model.mix)[2:NumVar], -sum(fixef(model.mix)[2:NumVar]))
result <- cbind(result, "LSMixed"=BLUEs)

#BLUPs NO. Non torna il numero!!!!!
#options(contrasts=c('contr.sum','contr.poly'))
#model.mix <- lme(Yield ~ e, random = ~ 1|v)
#BLUPs <- fixef(model.mix)[1] + ranef(model.mix)
#names(BLUPs) <- "BLUPs"
#print(as.vector(ranef(model.mix)))
#result <- cbind(result, BLUPs)
return(result)
}

YieldIndices <- function(Yield, Var, Env){
#Calcolo degli indici produttivi, su prove multiambiente
#Lavora sulle medie (tabelle a due vie genotipo x ambiente)  
#Data: 15/11/2009

MedieEnv <- tapply(Yield, Env, mean)[Env]
Index <- as.vector(Yield/MedieEnv*100)
Index
}

stability.met <- function(Yield, Var, Env){
# Questa funzione calcola le statistiche descrittive per prove varietali
# ripetute in molti ambienti, anche sbilanciate. ATTENZIONE: lavora sulle 
# medie ambientali delle varieta (non sulle repliche). calcola gli indici
# produttivi e di stabilita
#Data: 25/11/2009

MedieEnv <- tapply(Yield, Env, mean)[Env]
MedieVar <- tapply(Yield, Var, mean)[Var]
MediaGen <- mean(Yield)
NumVar <- length(tapply(Yield, Var, mean))
NumEnv <- length(tapply(Yield, Env, mean))

Index <- as.vector(Yield/MedieEnv*100)
Ecoval <- tapply((Yield - MedieEnv - MedieVar + MediaGen)^2, Var, sum)
#Ecoval <- Ecoval/(NumEnv-1)
EnvirVar <- tapply(Yield, Var, var)
#msge <- EnvirVar*(NumEnv-1)
msge <- sum((Yield - MedieEnv - MedieVar + MediaGen)^2)
sigma <- 1/((NumVar-1)*(NumVar-2)*(NumEnv-1))*(NumVar*(NumVar-1)*Ecoval - msge)
result <- summary.Met(Index, Var, Env)[,1:5]

IndiciMagg100 <- Index
for(i in 1:length(Index)){   
      if(Index[i]>100) IndiciMagg100[i] <- 1 else IndiciMagg100[i] <- 0
}    
IndiciAlti <- tapply(IndiciMagg100, Var, sum)[row.names(result)]
IndiceAff <- result[,2]-0.675*result[,5]

result <- cbind(result, "Num>100"=IndiciAlti)
result <- cbind(result, "Affidab."=IndiceAff)
result <- cbind(result, "Ecoval."=Ecoval)
result <- cbind(result, "EnvirVar"=EnvirVar)
result <- cbind(result, "Shukla's s2" = sigma)
result
}
