biplot.met <- function(obj, biplot = 1, xlim=NULL, ylim=NULL, elabels=NULL, 
                       glabels=NULL, quad=F, percDev=F, size=2, cexG=0.9, cexE=0.9,
                       xlab=NULL, ylab=NULL, font=1){
    if (class(obj) == "GGEobject") biplot <- 2                  
    E <- obj$Environment_scores
    G <- obj$Genotype_scores
    print(G)
    envir.mean <- obj$Environment_means
    var.mean <- obj$Genotype_means
    overall.mean <- mean(envir.mean)
    if (biplot == 1) {
        plot(1, type = "n", 
            xlim = if(is.null(xlim)) range(c(envir.mean, var.mean)) else xlim, 
            ylim = if(is.null(ylim)) range(c(E[, 1], G[, 1])) else ylim, 
            xlab = if(is.null(xlab)) "VarY" else xlab, 
            ylab = if(is.null(ylab)) "PC 1" else ylab,)
        points(envir.mean, E[, 1], "n", col = "red", lwd = 5)
        text(envir.mean, E[, 1],
            labels = if(is.null(elabels)) sub(' +$', '', as.vector(row.names(envir.mean))) else elabels, 
            #adj = c(0.5, 0.5), 
            col = "red", font=font, 
            cex = cexE)
        points(var.mean, G[, 1], "n", col = "blue", lwd = 5)
        text(var.mean, G[, 1], 
            labels = if(is.null(glabels)) sub(' +$', '', as.vector(row.names(var.mean))) else glabels, 
            #adj = c(0.5, 0.5), 
            col = "blue", font=font, 
            cex = cexG)
        abline(h = 0, lty = 5)
        abline(v = overall.mean, lty=5) 
    }
    else {
        plot(1, type = "n", 
        xlim = if(is.null(xlim)) range(c(E[, 1], G[, 1])) else xlim, 
        ylim = if(is.null(ylim)) range(c(E[, 2], G[, 2])) else ylim,
        xlab = "PC 1", ylab = "PC 2")
        points(E[, 1], E[, 2], "n", col = "red", lwd = 5)
        text(E[, 1], E[, 2], 
        	labels = if(is.null(elabels)) sub(' +$', '', as.vector(row.names(E))) else elabels,
        	col = "red", font=2, 
        	cex=cexE)
        points(G[, 1], G[, 2], "n", col = "blue", lwd = 5)
        text(G[, 1], G[, 2], 
        labels = if(is.null(glabels)) sub(' +$', '', as.vector(row.names(G))) else glabels, 
        	#adj = c(0.5, 0.5), 
        	col = "blue", font=font, 
        	cex=cexG)
        if (quad==TRUE){
        abline(h = 0, lty = 5)
        abline(v = 0, lty=5)
        }
        if (percDev != FALSE){
          expvar <- c(paste(round(obj$GGE_summary$Perc_of_Total_SS[1]),"%"), paste(round(obj$GGE_summary$Perc_of_Total_SS[2]),"%"))
          arrows(percDev[1], percDev[2], (percDev[1]+max(c(E[, 1], G[, 1]))/7), percDev[2], length=0.15)
          text((percDev[1]+max(c(E[, 1], G[, 1]))/7), percDev[2], label=expvar[1], pos=4)
          arrows(percDev[1], percDev[2],percDev[1], (percDev[2]+max(c(E[, 1], G[, 1]))/7), length=0.15)
          text(percDev[1], (percDev[2]+max(c(E[, 1], G[, 1]))/7), label=expvar[2], pos=3)
        }
    }
}

segment <- function(x1, y1, x2, y2, lty=NULL){
  curve(((x - x1)*(y2 - y1)) / (x2 - x1) + y1, from=x1, to=x2, lty = if(is.null(lty)) 1 else lty, add=T)
}

perp2p <- function(x1, y1, x2, y2, x3, y3, lty=NULL){
  m <- (y2 - y1)/(x2 - x1)
  curve(-1/m*(x-x3)+y3, from=if(max(x1,x2)>=x3 & min(y1,y2)<=y3) x3, to=if(min(x1,x2)<=x3 & max(y1,y2)>=y3) x3, lty = if(is.null(lty)) 1 else lty, add=T) 
}
biplot.polygon <- function(vertex, obj){
E <- obj$Environment_scores
G <- obj$Genotype_scores
  for (i in 1:(length(vertex)-1)){
  segment(G[vertex[i],1], G[vertex[i],2], G[vertex[i+1],1], G[vertex[i+1],2], lty=2)
  #perp2p(G[(vertex[i]),1], G[(vertex[i]),2], G[(vertex[i+1]),1], G[(vertex[i+1]),2], 0, 0)
}
  segment(G[vertex[1],1], G[vertex[1],2], G[vertex[length(vertex)],1], G[vertex[length(vertex)],2], lty=2)
  #perp2p(G[vertex[1],1], G[vertex[1],2], G[vertex[length(vertex)],1], G[vertex[length(vertex)],2], 0, 0)
}

summary.Met <- function(Yield, Var, Env){
#Questa funzione calcola le statistiche descrittive per prove varietali
#ripetute in più ambienti, anche sbilanciate. ATTENZIONE: lavora sulle 
#medie ambientali delle varietà (non sulle repliche).
#Data: 15/04/2009

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

Index <- function(Yield, Var, Env){
#Questa funzione calcola le statistiche descrittive per prove varietali
#ripetute in più ambienti, anche sbilanciate. ATTENZIONE: lavora sulle 
#medie ambientali delle varietà (non sulle repliche). calcola medie
# max min SD Least squares means and Blups
#Data: 15/11/2009

MedieEnv <- tapply(Yield, Env, mean)[Env]
Index <- as.vector(Yield/MedieEnv*100)
Index
}

stability.Met <- function(Yield, Var, Env){
#Questa funzione calcola le statistiche descrittive per prove varietali
#ripetute in più ambienti, anche sbilanciate. ATTENZIONE: lavora sulle 
#medie ambientali delle varietà (non sulle repliche). calcola gli indici
# produttivi e di stabilità
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
msge <- EnvirVar*(NumEnv-1)
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
