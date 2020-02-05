#Exponential decay ###############################################
expoDecay.fun <- function(predictor, C0, k) {
  C0 * exp(- k * predictor)
}

expoDecay.Init <- function(mCall, LHS, data) {
  xy <- sortedXyData(mCall[["predictor"]], LHS, data)
  lmFit <- lm(log(xy[, "y"]) ~ xy[, "x"])
  coefs <- coef(lmFit)
  C0 <- exp(coefs[1])
  k <- - coefs[2]
  value <- c(C0, k)
  names(value) <- mCall[c("C0", "k")]
  value
}

NLS.expoDecay <- selfStart(expoDecay.fun, expoDecay.Init, parameters=c("C0", "k"))

"DRC.expoDecay" <-
  function(fixed = c(NA, NA), names = c("init", "k"))
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
      
      init <- parmMat[, 1]; k <- parmMat[, 2]
      init * exp ( - k * x )
    }
    
    ## Defining self starter function
    ssfct <- function(dataf)
    {
      x <- dataf[, 1]
      y <- dataf[, 2]
      
      ## Linear regression on pseudo y values
      pseudoY <- log( y + 0.000001)
      coefs <- coef( lm(pseudoY ~ x) )
      init <- exp(coefs[1])
      m <- - coefs[2]
      
      return(c(init, m)[notFixed])
    }
    
    ## Defining names
    pnames <- names[notFixed]
    
    ## Defining derivatives
    
    ## Defining the ED function
    
    ## Defining the inverse function
    
    ## Defining descriptive text
    text <- "Exponential Decay Model"
    
    ## Returning the function with self starter and names
    returnList <- list(fct = fct, ssfct = ssfct, names = pnames, text = text, noParm = sum(is.na(fixed)))
    
    class(returnList) <- "drcMean"
    invisible(returnList)
  }
