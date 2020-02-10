#Logarithmic regression ########################################
logCurve.fun <- function(predictor, a, b) {
  x <- predictor
  a  + b * log(x)
}

logCurve.Init <- function(mCall, LHS, data) {
  xy <- sortedXyData(mCall[["predictor"]], LHS, data)
  x <-  xy[, "x"]; y <- xy[, "y"]
  pseudoY <-  y 
  pseudoX <- log(x)
  coefs <- coef( lm(pseudoY ~ pseudoX) )
  a <- coefs[1]        
  b <- coefs[2]
  value <- c(a, b)
  names(value) <- mCall[c("a", "b")]
  value
}

NLS.logCurve <- selfStart(logCurve.fun, logCurve.Init, parameters=c("a", "b"))

"DRC.logCurve" <-
  function(fixed = c(NA, NA), names = c("a", "b"))
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
      a  + b * log(x)
    }
    
    ## Defining self starter function
    ssfct <- function(dataf)
    {
      x <- dataf[, 1]
      y <- dataf[, 2]
      
      #regression on pseudo y values
      pseudoY <-  y
      pseudoX <- log(x)
      coefs <- coef( lm(pseudoY ~ pseudoX) )
      a <- coefs[1]
      b <- coefs[2]
      
      return(c(a, b)[notFixed])
    }
    
    ## Defining names
    pnames <- names[notFixed]
    
    ## Defining derivatives
    
    ## Defining the ED function
    
    ## Defining the inverse function
    
    ## Defining descriptive text
    text <- "Linear regression on log-transformed x"
    
    ## Returning the function with self starter and names
    returnList <- list(fct = fct, ssfct = ssfct, names = pnames, text = text, noParm = sum(is.na(fixed)))
    
    class(returnList) <- "drcMean"
    invisible(returnList)
  }

#Logarithmic regression without intercept ################################
logCurveNI.fun <- function(predictor, b) {
  x <- predictor
  b * log(x)
}

logCurveNI.Init <- function(mCall, LHS, data) {
  xy <- sortedXyData(mCall[["predictor"]], LHS, data)
  x <-  xy[, "x"]; y <- xy[, "y"]
  pseudoY <-  y 
  pseudoX <- log(x)
  coefs <- coef( lm(pseudoY ~ pseudoX - 1) )
  b <- coefs[1]
  value <- c(b)
  names(value) <- mCall[c("b")]
  value
}

NLS.logCurveNI <- selfStart(logCurveNI.fun, logCurveNI.Init, parameters=c("b"))
