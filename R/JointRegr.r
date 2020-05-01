JointRegr <- function(Block, Var, Loc, Yield){
    #Richiede che il disegno sia bilanciato!!!!!
    #Returns: Regression slope BETA (must be different from 0). Eberhart and Russel 1966
    #              Regression slope b (must be different from 1) Finlay
    #              S2di: deviation from regression (Eberhart and Russel, 1966)
    #              each indicator is followed by significance test and related probability
    
    #Example from Jawahar R. Sharma, 2006. 
    #Statistical And Biometrical Techniques In Plant Breeding
    #New Age International ltd. New Delhi, India
    #Block <- rep(c(1:3), 42)
    #Var <- rep(LETTERS[1:7],each=18)
    #Loc <- rep(rep(letters[1:6], each=3), 7)
    #Yield <- c(60,65,60,80,65,75,70,75,70,72,82,90,48,45,50,50,40,40,
    #           80,90,83,70,60,60,85,90,90,70,85,80,40,40,40,38,40,50,
    #           25,28,30,40,35,35,35,30,30,40,35,35,35,25,20,35,30,30,
    #           50,65,50,40,40,40,48,50,52,45,45,50,50,50,45,40,48,40,
    #           52,50,55,55,54,50,40,40,60,48,38,45,38,30,40,35,40,35,
    #           22,25,25,30,28,32,28,25,30,26,28,28,45,50,45,50,50,50,
    #           30,30,25,28,34,35,40,45,35,30,32,35,45,35,38,44,45,40)
     
    #ANOVA
    options(contrasts=c('contr.sum','contr.poly'))
    options(contrasts=c('contr.treatment','contr.poly'))
    #Loc <- datiFull$Loc; Block <- datiFull$Block
    #Var <- datiFull$Var; Yield <- datiFull$Yield
    model <- aov(Yield ~ factor(Loc)/factor(Block) + factor(Var) + factor(Loc) + factor(Var):factor(Loc))
    
    #Environment random - better do this fixed!!!!!
    aovTable <- summary(model)[[1]]
    aovTable[2, 4] <-  aovTable[2, 3]/aovTable[4, 3]
    aovTable[2, 5] <-  pf(aovTable[2, 4], aovTable[2, 1], aovTable[4, 1], lower.tail=FALSE)
    #print(aovTable)
    
    MedieLoc <- tapply(Yield, Loc, mean)
    MedieVar <- tapply(Yield, Var, mean)
    GrandMean <- mean(Yield)
    
    DataRegr <- data.frame(Yield, Block, Var, Loc)
    DataRegr <- DataRegr[order(Loc),]
    MedieLoc <- rep(MedieLoc, each=length(levels(factor(Block)))*length(levels(factor(Var))))
    
    DataRegr <- cbind(DataRegr, MedieLoc)
    DataRegr <- DataRegr[order(Var),]
    MedieVar <- rep(MedieVar, each=length(levels(factor(Block)))*length(levels(factor(Loc))))
    
    DataRegr <- cbind(DataRegr, MedieVar, EnvIndex=MedieLoc-GrandMean)
    DataRegr <- DataRegr[order(Loc),]
    DataRegr$Locf <- factor(DataRegr$Loc)
    
    Joint <- lm(Yield ~ factor(Var)/MedieLoc, data=DataRegr)
    aovJoint <- anova(Joint)
    
    #Perkins-Jinks
    EteroSS <- aovJoint[2,2] - aovTable[1,2]
    EteroDF <- aovJoint[2,1] - 1
    EteroMS <- EteroSS/EteroDF
    ResSS <- aovTable[4,2] - EteroSS
    ResDF <- aovTable[4,1] - EteroDF
    ResMS <- ResSS/ResDF
    row1 <- aovTable[2,]
    row2 <-  aovTable[1,]
    row3 <-  aovTable[4, ]
    row4 <- row3
    row.names(row4) <- "   Joint Regression"
    row4[1] <- EteroDF; row4[2] <- EteroSS; row4[3] <- EteroMS; row4[4] <-  EteroMS/ResMS;  row4[5] <-  pf(as.numeric(row4[4]), as.numeric(row4[1]), as.numeric(ResDF), lower.tail=FALSE)
    row5 <- row3
    row.names(row5) <- "   Deviation"
    row5[1] <- ResDF; row5[2] <- ResSS; row5[3] <- ResMS; row5[4] <- NA; row5[5] <- NA; 
    row6 <-  aovTable[5, ]
    row7 <- aovTable[3,]
    row.names(row7) <- "Blocks in Environments"
    aovJoint <- rbind(row2, row7, row1, row3, row4, row5, row6)
    PooledError <- row6[1,3]
    RSquare <- as.numeric(row4[2]/row3[2])
  
    
    #Calcolo coefficienti
    DataTemp <- aggregate(list(DataRegr$Yield, DataRegr$MedieLoc, DataRegr$MedieVar), by=list(DataRegr$Locf,DataRegr$Var), mean)
    names(DataTemp)[1] <- "Env"; names(DataTemp)[2] <- "Gen"; names(DataTemp)[3] <- "Yield"; names(DataTemp)[4] <- "LocMeans"; names(DataTemp)[5] <- "GenMeans"; 
    GEeff <- DataTemp$Yield - DataTemp$LocMeans - DataTemp$GenMeans + mean(DataTemp$Yield) 
    DataTemp <- data.frame(DataTemp, GEeff)
    #print(DataTemp[order(DataTemp$Env),])
    
    #i <- 2
    beta <- c(); b <- c(); SSreg <- c(); SSdev <- c(); MSdev <- c(); F <- c(); ProbF <- c(); Sd2 <- c(); FDev <- c(); ProbFDev <- c() 
    #print( length(tapply(Yield, Var, mean)))
    for(i in 1:length(tapply(Yield, Var, mean))){
        tempF <- subset(DataTemp, Gen==levels(factor(DataTemp$Gen))[i])
        prova <- mean(tempF$LocMeans)
        #temp <- lm(GEeff ~ I(LocMeans-prova), data=tempF)
        temp <- lm(Yield ~ I(LocMeans-prova), data=tempF)
        
        #print(anova(temp))
        beta[i] <- coef(temp)[2]
        SSreg[i] <- anova(temp)[1,2]
        SSdev[i] <- anova(temp)[2,2]
        MSdev[i] <- anova(temp)[2,3]
        F[i] <- anova(temp)[1,4]
        ProbF[i] <- anova(temp)[1,5]
        b[i] <- beta[i] - 1
        Sd2[i] <- MSdev[i]-PooledError/length(levels(factor(Block)))
        FDev[i] <- MSdev[i]/(PooledError/length(levels(factor(Block))))
        #print(anova(temp)[1,3]); print(row cfhj6[1,1])
        ProbFDev[i] <- 1 - pf(FDev[i], anova(temp)[2,1], row6[1,1])
    }

  JointTable <- cbind("MedieVar" = tapply(Yield, Var, mean), "BETA"=beta, "b" = b, "SSReg" = SSreg, "SSDev" = SSdev, "MSdev" = MSdev, "F" = F, "ProbF"=ProbF, "Sd2" = Sd2, "FDev" = FDev, "ProbFDev" = ProbFDev)
  finRes <- list("anova" = aovJoint, "DeterminationCoefficient" = RSquare, "JointAnalysis" = JointTable)
  return(finRes)
}


JointRegrM <- function(Var, Loc, Yield, sigma.e = 0){
  # Richiede che il disegno sia bilanciato!!!!!
  # Returns: Regression slope BETA (must be different from 0). Eberhart and Russel 1966
  # Regression slope b (must be different from 1) Finlay
  # S2di: deviation from regression (Eberhart and Russel, 1966)
  # need an estimate of the RSS from ANOVA on replicates
  # each indicator is followed by significance test and related probability
  
  Var <- factor(Var)
  Loc <- factor(Loc)
  #ANOVA
  library(nlme)
  CovL <- fitted(lm(Yield ~ Loc-1)) - mean(Yield)
  CovV <- fitted(lm(Yield ~ Var-1)) - mean(Yield)
  Joint <- gls(Yield ~ Var/CovL - 1, weights=varIdent(form=~1|Var))
  summary(Joint)$tTable
  aovJoint <- anova(Joint)
  
  GEeff <- Yield - fitted(lm(Yield ~ Loc-1)) - fitted(lm(Yield ~ Var-1)) + mean(Yield) 
  DataTemp <- data.frame(Yield, Var, Loc, CovL, CovV, GEeff)
  
  beta <- c(); b <- c(); SSreg <- c(); SSdev <- c(); MSdev <- c(); F <- c(); ProbF <- c(); Sd2 <- c(); FDev <- c(); ProbFDev <- c() 
  #print( length(tapply(Yield, Var, mean)))
  for(i in 1:length(tapply(Yield, Var, mean))){
      tempF <- subset(DataTemp, Var==levels(factor(DataTemp$Var))[i])
      prova <- mean(tempF$CovL)
      temp <- lm(Yield ~ I(CovL-prova), data=tempF)
      beta[i] <- coef(temp)[2]
      SSreg[i] <- anova(temp)[1,2]
      SSdev[i] <- anova(temp)[2,2]
      MSdev[i] <- anova(temp)[2,3]
      F[i] <- anova(temp)[1,4]
      ProbF[i] <- anova(temp)[1,5]
      b[i] <- beta[i] - 1
      Sd2[i] <- MSdev[i]
      }
  
  JointTable <- cbind("BETA"=beta, "b" = b, "SSReg" = SSreg, "SSDev" = SSdev, "MSdev" = MSdev, "F" = F, "ProbF"=ProbF, "Sd2" = Sd2 - sigma.e, "FDev" = FDev, "ProbFDev" = ProbFDev)
  row.names(JointTable) <- levels(Var)
  finRes <- list("GLS" = summary(Joint)$tTable, "JointAnalysis" = JointTable)
  return(finRes)
}
