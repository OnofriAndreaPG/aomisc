# Exponential growth function ##################
expoGrowth.fun <- function(predictor, init, k) {
  init * exp(k * predictor)
}

expoGrowth.Init <- function(mCall, LHS, data, ...) {
  xy <- sortedXyData(mCall[["predictor"]], LHS, data)
  x <-  xy[, "x"]; y <- xy[, "y"]
  ## Linear regression on pseudo y values
  pseudoY <- log(y + 0.000001)
  coefs <- coef( lm(pseudoY ~ x) )
  k <- coefs[1]
  b <- coefs[2]
  a <- exp(k)
  value <- c(a, b)
  names(value) <- mCall[c("init", "k")]
  value
}

NLS.expoGrowth <- selfStart(expoGrowth.fun, 
                            expoGrowth.Init, 
                            parameters=c("init", "k"))

"DRC.expoGrowth" <-
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
      init * exp ( k * x )
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
      m <- coefs[2]
      
      return(c(init, m)[notFixed])
    }
    
    ## Defining names
    pnames <- names[notFixed]
    
    ## Defining derivatives
    deriv1 <- function(x, parms){
      parmMat <- matrix(parmVec, nrow(parms), 
                        numParm, byrow = TRUE)
      parmMat[, notFixed] <- parms
      
      # Approximation by using finite differences
      a <- as.numeric(parmMat[,1])
      b <- as.numeric(parmMat[,2])
      
      d1.1 <- expoGrowth.fun(x, a, b)
      d1.2 <- expoGrowth.fun(x, (a + 10e-7), b)
      d1 <- (d1.2 - d1.1)/10e-7
      
      d2.1 <- expoGrowth.fun(x, a, b)
      d2.2 <- expoGrowth.fun(x, a, (b + 10e-7) )
      d2 <- (d2.2 - d2.1)/10e-7
      
      cbind(d1, d2)[notFixed]
    }
    
    ## Defining the first derivative (in x=dose)
    ##  based on deriv(~c+(d-c)*(exp(-exp(b*(log(x)-log(e))))), "x", function(x, b,c,d,e){})
    derivx <- function(x, parm)
    {
      parmMat <- matrix(parmVec, nrow(parm), numParm, byrow = TRUE)
      parmMat[, notFixed] <- parm
      
      a <- as.numeric(parmMat[,1])
      b <- as.numeric(parmMat[,2])
      
      d1.1 <- expoGrowth.fun(x, a, b)
      d1.2 <- expoGrowth.fun((x + 10e-7), a, b)
      d1 <- (d1.2 - d1.1)/10e-7
      
      d1
      
    }
    
    ## Defining the ED function
    
    ## Defining the inverse function
    
    ## Defining descriptive text
    text <- "Exponential Growth Model"
    
    ## Returning the function with self starter and names
    returnList <- list(fct = fct, ssfct = ssfct, 
                       deriv1 = deriv1, derivx = derivx,
                       names = pnames, text = text, 
                       noParm = sum(is.na(fixed)))
    
    class(returnList) <- "drcMean"
    invisible(returnList)
  }


