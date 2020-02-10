#Yield-Loss function (rectangular hyperbola) ######################
YL.fun <- function(predictor, i, A) {
  i * predictor/(1 + i/A * predictor)
}

YL.Init <- function(mCall, LHS, data) {
  xy <- sortedXyData(mCall[["predictor"]], LHS, data)
  x <- xy[, "x"]; y <- xy[, "y"]
  pseudoX <- 1 / x[x > 0]; pseudoY <- 1 / y[x > 0]
  lmFit <- lm(pseudoY ~ pseudoX)
  coefs <- coef(lmFit)
  A <- 1 / coefs[1]
  i <- 1 / coefs[2]
  value <- c(i, A)
  names(value) <- mCall[c("i", "A")]
  value
}
NLS.YL <- selfStart(YL.fun, YL.Init, parameters=c("i", "A"))

"DRC.YL" <- function(fixed = c(NA, NA), names = c("i", "A")) {
  ## Checking arguments
  numParm <- 2
  if (!is.character(names) | !(length(names) == numParm)) {stop("Not correct 'names' argument")}
  if (!(length(fixed) == numParm)) {stop("Not correct 'fixed' argument")}
  
  ## Fixing parameters (using argument 'fixed')
  notFixed <- is.na(fixed)
  parmVec <- rep(0, numParm)
  parmVec[!notFixed] <- fixed[!notFixed]
  
  ## Defining the non-linear function
  fct <- function(x, parm) {
    parmMat <- matrix(parmVec, nrow(parm), numParm, byrow = TRUE)
    parmMat[, notFixed] <- parm
    i <- parmMat[, 1]; A <- parmMat[, 2]
    YL.fun(x, i, A)
  }
  
  ## Defining self starter function
  ssfct <- function(dataf) {
    x <- dataf[, 1]
    y <- dataf[, 2]
    
    #regression on pseudo y values
    pseudoY <- 1 /  y[x > 0] 
    pseudoX <- 1 / x [x > 0] 
    coefs <- coef( lm(pseudoY ~ pseudoX) )
    A <- 1 / coefs[1]; i <- 1 / coefs[2]
    
    return(c(i, A)[notFixed])
  }
  
  ## Defining names
  pnames <- names[notFixed]
  
  ## Defining derivatives
  
  ## Defining the ED function
  
  ## Defining the inverse function
  
  ## Defining descriptive text
  text <- "Yield-Loss function (Cousens, 1985)"
  
  ## Returning the function with self starter and names
  returnList <- list(fct = fct, ssfct = ssfct, names = pnames, text = text, noParm = sum(is.na(fixed)))
  
  class(returnList) <- "drcMean"
  invisible(returnList)
}
