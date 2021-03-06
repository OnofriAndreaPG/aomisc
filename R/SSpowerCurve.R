#Power Curve ########################################################
powerCurve.fun <- function(predictor, a, b) {
  a * ( predictor ^ b )
}

powerCurve.Init <- function(mCall, LHS, data) {
  xy <- sortedXyData(mCall[["predictor"]], LHS, data)
  lmFit <- lm(log(xy[, "y"]) ~ log(xy[, "x"]))
  coefs <- coef(lmFit)
  a <- exp(coefs[1])
  b <- coefs[2]
  value <- c(a, b)
  names(value) <- mCall[c("a", "b")]
  value
}

NLS.powerCurve <- selfStart(powerCurve.fun, powerCurve.Init, parameters=c("a", "b"))

"DRC.powerCurve" <- function(fixed = c(NA, NA), names = c("a", "b"))
{
  ## Checking arguments
  numParm <- 2
  if (!is.character(names) | !(length(names) == numParm)) {stop("Not correct 'names' argument")}
  if (!(length(fixed) == numParm)) {stop("Not correct 'fixed' argument")}
  
  ## Fixing parameters (using argument 'fixed')
  notFixed <- is.na(fixed)
  parmVec <- rep(0, numParm)
  parmVec[!notFixed] <- fixed[!notFixed]
  
  ## Defining the non-linear function
  fct <- function(x, parm)
  {
    
    
    parmMat <- matrix(parmVec, nrow(parm), numParm, byrow = TRUE)
    parmMat[, notFixed] <- parm
    
    a <- parmMat[, 1]; b <- parmMat[, 2]
    a * x ^(b)
  }
  
  ## Defining self starter function
  ssfct <- function(dataf)
  {
    x <- dataf[, 1]
    y <- dataf[, 2]
    
    #regression on pseudo y values
    pseudoY <- log( y + 0.00001)
    pseudoX <- log(x)
    coefs <- coef( lm(pseudoY ~ pseudoX) )
    a <- exp(coefs[1])
    
    b <- coefs[2]
    
    return(c(a, b)[notFixed])
  }
  
  ## Defining names
  pnames <- names[notFixed]
  
  ## Defining derivatives
  
  ## Defining the ED function
  
  ## Defining the inverse function
  
  ## Defining descriptive text
  text <- "Power curve (Freundlich equation)"
  
  ## Returning the function with self starter and names
  returnList <- list(fct = fct, ssfct = ssfct, names = pnames, text = text, noParm = sum(is.na(fixed)))
  
  class(returnList) <- "drcMean"
  invisible(returnList)
}

