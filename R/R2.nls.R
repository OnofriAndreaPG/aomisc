R2nls <- function(object){
  
  if (!inherits(object, "nls") & !inherits(object, "drc"))
    stop("use only with \"nls\" or \"drc\" objects")
  
  if (inherits(object, "nls")){
    formula <- as.formula(summary(object)$formula)
    varNames <- all.vars(formula)
    if(!is.environment(eval(object$data))){
      Y <- eval(object$data)[,varNames[1]]
    } else {
      Y <- get(varNames[1])
    }
    dfMod <- summary(object)$df[2]
    sigma <- summary(object)$sigma
  } else if (inherits(object, "drc")){
    # object$origData
    varNames <- all.vars(object$call$formula)
    Y <- eval(object$data)[,varNames[1]]
    dfMod <- summary(object)$df
    sigma <- sqrt( summary(object)$resVar )
  }
  Y <- as.numeric(unlist(Y))
  SSm <- sum( (fitted(object) - mean(Y))^2) # Regression deviance
  SSt <- deviance( lm(Y ~ 1) ) # Total deviance about the mean
  SSr <- deviance(object)  # Residual deviance
  PseudoR2 <- (SSt - SSr)/SSt # R2 as in Schebenberger Eq. 5.24 (pag. 211)
  R2 <- SSm/SSt # R2 traditional as in Scebenberger Eq. 5.23
  
  # MSt <- SSt/(length(Y) - 1)
  # MSr <- SSr/dfMod
  # R2adj <- 1 - MSr/MSt # Adjusted R2
  # 
  # MSE <- SSr/dfMod
  # RMSE <- sigma
  # RRMSE <- RMSE/mean(Y)
  
  # R2 generalised: to be done for GLMs
  # ll1 <- as.numeric(logLik(object))
  # ll2 <- as.numeric(logLik( lm(Y ~ 1) ))
  # R2gen <- 1 - exp(-2/length(Y)*(ll1 - ll2)) # Schabenberger, 6.46 pag 343
  # R2genResc <- R2gen/(1 - exp(2/length(Y) * ll2)) # Schabenberger, 6.48 pag 344
  returnList <- list(PseudoR2 = PseudoR2, R2 = R2 
                     # R2adj = R2adj,
                     #R2gen = R2gen, R2gen.rescaled = R2genResc,
                     # , MSE = MSE, RMSE = RMSE, RRMSE = RRMSE
                     )
  returnList
  }
