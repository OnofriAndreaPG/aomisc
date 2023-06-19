## GENOTYPE PLUS GENOTYPE X ENVIRONMENT INTERACTION, 
## AS APPLIED TO AGRICULTURE VARIETY TRIALS
## REFERENCE: Yan et al., 2000. Cultivar evaluation and mega-environment investigation
## based on the GGE biplot
## Update: 26/05/2023

## INPUTS:
## variety is a vector of variety codes (factor)
## envir is a vector of environment codes (factor)
## block is a vector of block codes (factor)
## yield is a vector of yields to be analysed (numerical)
## PC is the number of PC to be considered (default is 2)
## GGE: works on balanced raw data
## GGEmeans: works on means
## scale is the partitioning of singular values

GGE <- function(yield, genotype, environment, block = NULL,
                 PC = 2) {
  MSE <- NULL
  dfr <- NULL

## 1 - Model fitting and descriptive statistics
  
  envir <- as.factor(environment)
  variety <- as.factor(genotype)
  
  # Verify whether it is balanced and replicated
  conTab <- table(envir, variety)
  if(any(conTab == 0)) stop("Missing cells were found. Imputing values is necessary")
  balanced <- ifelse(length(unique(conTab)) == 1, TRUE, FALSE)
  unreplicated <- FALSE
  if(balanced) unreplicated <- ifelse(unique(conTab) == 1, TRUE, FALSE)
  numReps <- as.numeric(unique(conTab))
  if(!is.null(MSE) & !is.null(dfr)) externalError <- TRUE else externalError <- FALSE
  
  # Estimate the MSE if not passed in (not sure it is useful) ...
  add.anova <- NULL
  
  if(!is.null(block) & is.null(MSE) & !unreplicated) {
    # With blocks, it includes the block %in% environments in the model
    add.anova <- lm(yield ~ variety*envir + block %in% envir)
    suppressMessages( tab <- as.data.frame(emmeans::emmeans(add.anova, ~variety:envir)) )
    MSE <- summary(add.anova)$sigma^2
    MSE <- MSE/numReps
    df.res <- add.anova$df.residual
  } else if(is.null(block) & is.null(MSE) & !unreplicated) {
    # It does not include the blocks
    add.anova <- lm(yield ~ variety*envir)
    suppressMessages( tab <- as.data.frame(emmeans::emmeans(add.anova, ~variety:envir)) )
    # tab <- data.frame(variety, envir, emmean = yield)
    MSE <- summary(add.anova)$sigma^2
    MSE <- MSE/numReps
    df.res <- add.anova$df.residual
  } else if(is.null(MSE) & unreplicated) {
    add.anova <- lm(yield ~ variety*envir)
    suppressWarnings(suppressMessages( tab <- as.data.frame(emmeans::emmeans(add.anova, ~variety:envir)) ))
    MSE <- 0
    df.res <- 0
    
  } else {
    add.anova <- lm(yield ~ variety*envir)
    suppressWarnings(suppressMessages( tab <- as.data.frame(emmeans::emmeans(add.anova, ~variety:envir)) ))
    df.res <- add.anova$df.residual
  }
  
  if(is.null(dfr)) dfr <- df.res
  
  # Calculate GGE matrix (ge.eff) and means for main effects
  int.mean <- unstack(tab, form = emmean ~ envir)
  rownames(int.mean) <- unique(tab$variety)
  
  ge.eff <- as.data.frame(scale(int.mean, center=T, scale=F))
  
  suppressWarnings(suppressMessages({ 
    overall.mean <- as.data.frame(emmeans(add.anova, ~1))$emmean
    envir.mean <- as.data.frame(emmeans::emmeans(add.anova, ~envir))$emmean
    var.mean <- as.data.frame(emmeans::emmeans(add.anova, ~variety))$emmean
  }))
  names(var.mean) <- rownames(int.mean)
  names(envir.mean) <- colnames(int.mean)
  envir.num <- length(envir.mean)
  var.num <- length(var.mean)
  
  
## 2 - Variance partitioning (additive model)
if(!unreplicated) tab <- anova(add.anova) else tab <- NULL
  
## 3 - SVD
  scale <- 0.5
  dec <- svd(ge.eff, nu = PC, nv = PC)
  if (PC > 1){ 
    D <- diag(dec$d[1:PC]^scale)
    Dbis <- diag(dec$d[1:PC]^(1-scale)) 
  } else {
    D <- dec$d[1:PC]^scale
    Dbis <- dec$d[1:PC]^(1-scale)
  }
  G<-dec$u %*% D
  E<-dec$v %*% Dbis
  
  Ecolnumb <- c(1:PC)
  Ecolnames <- paste("PC", Ecolnumb, sep = "")
  dimnames(E) <- list(levels(envir), Ecolnames)
  dimnames(G) <- list(levels(variety), Ecolnames)
  
## 4 - Singular values and significance of PCs  
  svalues <- dec$d
  PC.n <- c(1:length(svalues))
  PC.SS <- svalues^2
  percSS <- PC.SS/sum(PC.SS)*100
  GGE.table <- data.frame("PC"=PC.n,"Singular_value"=svalues,"PC_SS"=PC.SS, "Perc_of_Total_SS"=percSS)
  numblock <-length(levels(block))
  # GGE.SS <- (t(as.vector(ge.eff)) %*% as.vector(ge.eff)) * numblock
  GGE.SS <- sum(ge.eff^2) * numblock
  
  # int.SS <- sum(int.eff^2)
  # PC.SS <-  (dec$d[1:PC]^2) # * numblock 
  # PC.DF <- var.num + envir.num - 1 - 2*Ecolnumb
  # residual.SS <- int.SS - sum(PC.SS)
  # residual.DF <- ((var.num - 1)*(envir.num - 1)) - sum(PC.DF)
  # PC.SS[PC + 1] <- residual.SS  
  # PC.DF[PC + 1] <- residual.DF
  # MS <- PC.SS/PC.DF
  
  
  # if(balanced | externalError){
  #   Fcomp <- MS/MSE
  #   probab <- pf(Fcomp, PC.DF, dfr, lower.tail = FALSE)
  #   percSS <- PC.SS/int.SS
  #   rowlab <- c(Ecolnames, "Residuals")
  #   mult.anova <- data.frame(Effect = rowlab, SS = PC.SS, DF = PC.DF, MS = MS, F = Fcomp, Prob. = probab, "Perc_of_Total_SS"=percSS)
  # } else {
  #   percSS <- PC.SS/int.SS
  #   rowlab <- c(Ecolnames, "Residuals")
  #   mult.anova <- data.frame(Effect = rowlab, SS = PC.SS, DF = PC.DF, MS = MS, "Perc_of_Total_SS"=percSS)
  #   MSE <- NULL
  #   dfr <- NULL
  # }
  # 
  # if(PC >=2) stability <- sqrt((PC.SS[1]/PC.SS[2]*G[,1])^2+G[,2]^2) else stability = NA
  
## 6 - Results  
  result <- list(genotype_means=var.mean, environment_means=envir.mean, interaction_means=int.mean, 
      GE_effect=ge.eff, additive_ANOVA = tab, "GGE_Sum of Squares"=GGE.SS, 
      GGE_summary = GGE.table, environment_scores=E, "genotype_scores"=G)
  class(result) <- "GGEobject"
     
  return(invisible(result))
}

