# Hill function
hillCurveMean <- function(predictor, a, b, c) {
  (a * predictor^c)/(b + predictor^c)
}

hillCurveInit <- function(mCall, LHS, data, ...) {
  xy <- sortedXyData(mCall[["predictor"]], LHS, data)
  x <-  xy[, "x"]; y <- xy[, "y"]
  a <- max(y) * 1.05
  pseudoY <-  log(( a - y )/ y)
  pseudoX <- log(x)
  lmFit <- lm(pseudoY ~ pseudoX )
  coefs <- coef(lmFit)
  b <- exp(coefs[1])
  c <- - coefs[2]
  value <- c(a, b, c)
  names(value) <- mCall[c("a", "b", "c")]
  value
}

NLShillCurve <- selfStart(hillCurveMean, hillCurveInit, parameters=c("a", "b", "c"))

"hill" <-
  function(fixed = c(NA, NA, NA), names = c("a", "b", "c"))
  {
    ## Checking arguments
    numParm <- 3
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
      
      a <- parmMat[, 1]; b <- parmMat[, 2]; c <- parmMat[, 3]
      (a * x^c)/(b + x^c)
    }
    
    ## Defining self starter function
    ssfct <- function(dataf)
    {
      x <- dataf[, 1]
      y <- dataf[, 2]
      
      a <- max(y) * 1.05
      
      ## Linear regression on pseudo y values
      pseudoY <-  log(( a - y )/ y)
      pseudoX <- log(x)
      coefs <- coef( lm(pseudoY ~ pseudoX ) )
      b <- exp(coefs[1])
      c <- - coefs[2]
      
      return(c(a, b, c)[notFixed])
    }
    
    ## Defining names
    pnames <- names[notFixed]
    
    ## Defining derivatives
    
    ## Defining the ED function
    
    ## Defining the inverse function
    
    ## Defining descriptive text
    text <- "Hill function (Morgan-Mercer-Flodin)"
    
    ## Returning the function with self starter and names
    returnList <- list(fct = fct, ssfct = ssfct, names = pnames, text = text, noParm = sum(is.na(fixed)))
    
    class(returnList) <- "drcMean"
    invisible(returnList)
  }

"hillMax" <-
  function(fixed = c(NA, NA, NA), names = c("a", "b", "c"))
  {
    ## Checking arguments
    numParm <- 3
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
      
      a <- parmMat[, 1]; b <- parmMat[, 2]; c <- parmMat[, 3]
      (8 * 0.8 * (1 + b * exp ( -c * log(a))))/(1 + b * exp ( -c * log(x)))
    }
    
    ## Defining self starter function
    
    ## Defining names
    pnames <- names[notFixed]
    
    ## Defining descriptive text
    text <- "Hill function (Morgan-Mercer-Flodin) with expected value parameters replacement"
    
    ## Returning the function with self starter and names
    returnList <- list(fct = fct, names = pnames, text = text, noParm = sum(is.na(fixed)))
    
    class(returnList) <- "drcMean"
    invisible(returnList)
  }

