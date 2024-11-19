AnscombeTukey <- function(mod, premium = 0.025){
  # Implementing the rejection rule of Anscombe and Tukey, 1963
  # It works iteratively
  # Rejection rule function
  AnscombeTukeyRule <- function(mod, premium){
  
  # It must be a linear model object
  if(!any(inherits(mod, "lm"), TRUE)) stop ("This function only works on linear model objects")

  ## model is a lm object
  epsilon <- mod$residuals
  ENNE <- length(epsilon)
  NI <- mod$df.residual
  prob <- (premium * NI / ENNE)
  zed <- abs(qnorm(prob))
  K <- 1.4 + 0.85 * zed
  CI <- K * (1 - (K^2 - 2) / (4 * NI)) * sqrt(NI / ENNE)
  devst <- summary(mod)$sigma
  maxteor <- CI * devst

  ## Find the position of an outlier
  MaxRes <- max(abs(epsilon))
  Pos <- as.numeric(which(abs(epsilon) >= maxteor))
  return(list(Pos, maxteor))
}
  positions <- c()
  i <- 1
  resList <- AnscombeTukeyRule(mod, premium)
  if(length(resList[[1]]) > 0) {
    rule <- TRUE
    newobj <- mod
    Pos <- resList[[1]]
    positions[i] <- Pos
    i <- i + 1
    } else rule <- FALSE 
  
  while(rule){
    newdata <- eval(newobj$call$data)[-Pos,]
    newobj <- update(newobj, data = newdata)
    resList <- AnscombeTukeyRule(newobj, premium)
    # print(resList)
    if(length(resList[[1]]) > 0) {
      rule <- TRUE
      # newobj <- mod
      Pos <- resList[[1]]
      positions[i] <- Pos
      i <- i + 1
    } else rule <- FALSE
  }
  return(positions)
}

