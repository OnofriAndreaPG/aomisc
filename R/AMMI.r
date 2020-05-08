## ADDITIVE MAIN EFFECTS MULTIPLICATIVE INTERACTION MODEL, AS APPLIED TO AGRICULTURE VARIETY TRIALS
## REFERENCE: H. F. Gollob, 1968. A statistical model which combine features of factor analytic and
## analysis of variance techniques. Psychometrika, 33, 1, 73-114.
## Update: 16/07/2019  

## INPUTS:
## genotype is a vector of genotype codes (factor)
## env is a vector of environment codes (factor)
## block is a vector of block codes (factor)
## yield is a vector of yields to be analysed (numerical)
## PC is the number of PC to be considered (default is 2)
## AMMI: works on balanced raw data
## AMMImeans: works on means
  
AMMImeans <- function(yield, genotype, environment, PC = 2,
                      MSE = NULL, dfr = NULL) {

  variety <- genotype; envir <- environment
  add.anova <- aov(yield ~ envir * variety)
  int.eff <- model.tables(add.anova, type = "effects", cterms = "envir:variety")$tables$"envir:variety"
  int.mean <- model.tables(add.anova, type = "means", cterms = "envir:variety")$tables$"envir:variety"
  envir.mean <- model.tables(add.anova, type = "means", cterms = "envir")$tables$"envir"
  var.mean <- model.tables(add.anova, type = "means", cterms = "variety")$tables$"variety"
  overall.mean <- model.tables(add.anova, type = "means")$tables[[1]]
  
## 3 - SVD
  dec <- svd(int.eff, nu = PC, nv = PC)
  if (PC > 1){ 
    D <- diag(dec$d[1:PC]) 
  } else {
    D <- matrix(dec$d[1:PC],1,1)
  }
  E <- dec$u %*% sqrt(D)
  G <- dec$v %*% sqrt(D)
  Ecolnumb <- c(1:PC)
  Ecolnames <- paste("PC", Ecolnumb, sep = "")
  dimnames(E) <- list(levels(envir), Ecolnames)
  dimnames(G) <- list(levels(variety), Ecolnames)
  #stability <- sqrt(G[,1]^2+G[,2]^2)

## 4 - Significance of PCs
  svalues<-dec$d
  PC.n<-c(1:length(svalues))
  PC.SS<-svalues^2
  percSS<-PC.SS/sum(PC.SS)*100
  GGE.table<-data.frame("PC"=PC.n,"Singular_value"=svalues,"PC_SS"=PC.SS, "Perc_of_Total_SS"=percSS,
                        cum_perc = cumsum(percSS))
  GGE.SS <- (t(as.vector(int.eff))%*%as.vector(int.eff))
    
  
  if(!is.null(MSE) & !is.null(dfr)){
    ngen <- length(int.mean[1,])
    nenv <- length(int.mean[,1])
    df_n <- ngen + nenv - 1 - 2 * GGE.table[,1]
    ss <- GGE.table$PC_SS
    ms <- ss/df_n
    f <- ms/MSE
    pf <- pf(f, df_n, dfr, lower.tail = F)
    aov.tab <- data.frame(PC = GGE.table[,1], SS = ss, DF = df_n, MS = ms,
                          "F" = f, "P value" = pf)
  } else {
    aov.tab <- NA
  }
    
## 6 - Other results  
  result <- list(means_table = int.mean,
       interaction_effect=int.eff, "AMMI_SS"=GGE.SS, 
       summary = GGE.table, anova = aov.tab,
       environment_scores = E, genotype_scores = G,
       "genotype_means"=var.mean, "environment_means"=envir.mean)#,Stability = stability)
  # cat(paste("Result of AMMI Analysis", "\n", "\n"))
  # class(result) <- "AMMIobject"
  # print(E)
  # cat(paste("\n","Genotype Scores", "\n"))
  # print(G)
  return(invisible(result))
}

AMMI <- function(yield, genotype, environment, block, PC = 2) {

## 1 - Descriptive statistics
  variety <- genotype; envir <- environment
  overall.mean <- mean(yield)
  envir.mean <- tapply(yield, envir, mean)
  var.mean <- tapply(yield, variety, mean)  
  int.mean <- tapply(yield, list(variety,envir), mean)
  envir.num <- length(envir.mean)
  var.num <- length(var.mean)
  
## 2 - ANOVA (additive model) 
  variety <- factor(variety)
  envir <- factor(envir)
  block <- factor(block)  
  add.anova <- aov(yield ~ envir + block %in% envir + variety + envir:variety)
  modelTables <- model.tables(add.anova, type = "effects", cterms = "envir:variety")
  int.eff <- modelTables$tables$"envir:variety"
  add.anova.residual.SS <- deviance(add.anova)
  add.anova.residual.DF <- add.anova$df.residual
  add.anova.residual.MS <- add.anova.residual.SS/add.anova.residual.DF
  anova.table <- summary(add.anova)
  row.names(anova.table[[1]]) <- c("Environments", "Genotypes", "Blocks(Environments)","Environments x Genotypes", "Residuals")
  
## 3 - SVD
  dec <- svd(int.eff, nu = PC, nv = PC)
  if (PC > 1){ 
    D <- diag(dec$d[1:PC]) 
  } else {
    D <- dec$d[1:PC]
  }
  E <- dec$u %*% sqrt(D)
  G <- dec$v %*% sqrt(D)
  Ecolnumb <- c(1:PC)
  Ecolnames <- paste("PC", Ecolnumb, sep = "")
  dimnames(E) <- list(levels(envir), Ecolnames)
  dimnames(G) <- list(levels(variety), Ecolnames)
  
## 4 - Significance of PCs
  numblock <- length(levels(block))
  int.SS <- c( (t(as.vector(int.eff)) %*% as.vector(int.eff))*numblock )
  PC.SS <-  (dec$d[1:PC]^2)*numblock  
  PC.DF <- var.num + envir.num - 1 - 2*Ecolnumb
  residual.SS <- int.SS - sum(PC.SS)
  residual.DF <- ((var.num - 1)*(envir.num - 1)) - sum(PC.DF)
  PC.SS[PC + 1] <- residual.SS  
  PC.DF[PC + 1] <- residual.DF
  MS <- PC.SS/PC.DF
  F <- MS/add.anova.residual.MS
  probab <- pf(F, PC.DF, add.anova.residual.DF, lower.tail = FALSE)
  percSS <- PC.SS/int.SS
  rowlab <- c(Ecolnames, "Residuals")
  mult.anova <- data.frame(Effect = rowlab, SS = PC.SS, DF = PC.DF, MS = MS, F = F, Prob. = probab, "% of Total SS"=percSS)
  stability <- sqrt(G[,1]^2+G[,2]^2)

## 6 - Results  
  result <- list(genotype_means = var.mean, environment_means = envir.mean,  interaction_means = int.mean,
       interaction_effect=int.eff, additive_ANOVA = anova.table, mult_Interaction = mult.anova, 
       environment_scores = E, genotype_scores = G, stability = stability)
       class(result) <- "AMMIobject"    
# cat(paste("Result of AMMI Analysis", "\n", "\n"))
# cat(paste("Environment Scores", "\n"))
# print(E)
# cat(paste("\n","Genotype Scores", "\n"))
# print(G)
return(invisible(result))
}

AMMI_old <- function(variety, envir, block, yield, PC=2, biplot=1) {

## 1 - Descriptive statistics
  overall.mean <- mean(yield)
  envir.mean <- tapply(yield, envir, mean)
  var.mean <- tapply(yield, variety, mean)  
  int.mean <- tapply(yield, list(variety,envir), mean)
  envir.num <- length(envir.mean)
  var.num <- length(var.mean)
  
## 2 - ANOVA (additive model) 
  variety <- factor(variety)
  envir <- factor(envir)
  block <- factor(block)  
  add.anova <- aov(yield ~ envir + block %in% envir + variety + envir:variety)
  modelTables <- model.tables(add.anova, type = "effects", cterms = "envir:variety")
  int.eff <- modelTables$tables$"envir:variety"
  add.anova.residual.SS <- deviance(add.anova)
  add.anova.residual.DF <- add.anova$df.residual
  add.anova.residual.MS <- add.anova.residual.SS/add.anova.residual.DF
  anova.table <- summary(add.anova)
  row.names(anova.table[[1]]) <- c("Environments", "Genotypes", "Blocks(Environments)","Environments x Genotypes", "Residuals")
  
## 3 - SVD
  dec <- svd(int.eff, nu = PC, nv = PC)
  if (PC > 1){ 
    D <- diag(dec$d[1:PC]) 
  } else {
    D <- dec$d[1:PC]
  }
  E <- dec$u %*% sqrt(D)
  G <- dec$v %*% sqrt(D)
  Ecolnumb <- c(1:PC)
  Ecolnames <- paste("PC", Ecolnumb, sep = "")
  dimnames(E) <- list(levels(envir), Ecolnames)
  dimnames(G) <- list(levels(variety), Ecolnames)
  
## 4 - Significance of PCs
  numblock <- length(levels(block))
  int.SS <- (t(as.vector(int.eff)) %*% as.vector(int.eff))*numblock
  PC.SS <- (dec$d[1:PC]^2)*numblock  
  PC.DF <- var.num + envir.num - 1 - 2*Ecolnumb
  residual.SS <- int.SS - sum(PC.SS)
  residual.DF <- ((var.num - 1)*(envir.num - 1)) - sum(PC.DF)
  PC.SS[PC + 1] <- residual.SS  
  PC.DF[PC + 1] <- residual.DF
  MS <- PC.SS/PC.DF
  F <- MS/add.anova.residual.MS
  probab <- pf(F, PC.DF, add.anova.residual.DF, lower.tail = FALSE)
  percSS <- PC.SS/int.SS
  rowlab <- c(Ecolnames, "Residuals")
  mult.anova <- data.frame(Effect = rowlab, SS = PC.SS, DF = PC.DF, MS = MS, F = F, Prob. = probab, "% of Total SS"=percSS)
  stability <- sqrt(G[,1]^2+G[,2]^2)

## 5 - Biplots
  if (biplot == 1) {
    plot(1, type = 'n', xlim = range(c(envir.mean, var.mean)), ylim = range(c(E[,1], G[,1])), xlab = "Yield",
    		 ylab = "PC 1")
    points(envir.mean, E[,1], "n", col = "red", lwd = 5)
    text(envir.mean, E[,1], labels = row.names(envir.mean), adj = c(0.5, 0.5), col = "red")
    points(var.mean, G[,1], "n", col = "blue", lwd = 5)
    text(var.mean, G[,1], labels = row.names(var.mean), adj = c(0.5, 0.5), col = "blue")
    abline(h = 0, v = overall.mean, lty = 5)
  } else {
    plot(1, type = 'n', xlim = range(c(E[,1], G[,1])), ylim = range(c(E[,2], G[,2])), xlab = "PC 1", 
    			ylab = "PC 2")
    points(E[,1], E[,2], "n",col = "red", lwd = 5)
    text(E[,1], E[,2], labels = row.names(E),adj = c(0.5,0.5),col = "red")
    points(G[,1],G[,2], "n", col = "blue", lwd = 5)
    text(G[,1], G[,2], labels = row.names(G),adj = c(0.5, 0.5), col = "blue")
  }
    
## 6 - Other results  
  result <- list(Genotype_means = var.mean, Environment_means = envir.mean,  Interaction_means = int.mean,
       Interaction_effect=int.eff, Additive_ANOVA = anova.table, Mult_Interaction = mult.anova, 
       Environment_scores = E, Genotype_scores = G, Stability = stability)
       class(result) <- "AMMIobject"
  return(result)
}

