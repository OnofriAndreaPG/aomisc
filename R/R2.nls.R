R2nls <- function(object){
  
  if (!inherits(object, "nls"))
    stop("use only with \"nls\" objects")

  formula <- as.formula(summary(object)$formula)
  varNames <- all.vars(formula)
  Y <- eval(object$data)[,varNames[1]]
  SSt <- deviance( lm(Y ~ 1) ) 
  SSr <- deviance(object)  
  R2 <- (SSt - SSr)/SSt
  MSt <- SSt/(length(Y) - 1)
  MSr <- SSr/summary(object)$df[2]
  R2adj <- 1 - MSr/MSt
  #1 - (1 - R2) * (length(Y) - 1) /(length(Y) - length(coef(object)) - 1)
  returnList <- list(R2 = R2, R2adj = R2adj)
  returnList
  }
