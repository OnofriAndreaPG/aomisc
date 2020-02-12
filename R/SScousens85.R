#Yield - Weed Density function ###################
cousens85.fun <- function(Dens, Ywf, i, A) {
  Ywf*(1 - (i * Dens) / (100 * (1 + i * Dens/A)))
}

"DRC.cousens85" <-
  function(fixed = c(NA, NA, NA), names = c("YWF", "i", "a"))
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
      
      YWF <- parmMat[, 1]; i <- parmMat[, 2]; a <- parmMat[, 3]
      YWF*(1 - (i*x) / (100 * (1 + i * x/a)))
    }
    
    ## Defining self starter function
    ssfct <- function(dataf)
    {
      x <- dataf[, 1]
      y <- dataf[, 2]
      
      YWF <- max(y)+10e-06
      YL <- (1 - y/YWF)*100
      #regression on pseudo y values
      pseudoY <- 1 /  (YL + 0.000001)
      pseudoX <- 1 / (x + 0.00001)
      coefs <- coef( lm(pseudoY ~ pseudoX) )
      a <- 1 / coefs[1]; i <- 1 / coefs[2]
      
      return(c(YWF, i, a)[notFixed])
    }
    
    ## Defining names
    pnames <- names[notFixed]
    
    ## Defining derivatives
    
    ## Defining the ED function
    
    ## Defining the inverse function
    
    ## Defining descriptive text
    text <- "Yield-Weed Density function (Cousens, 1985)"
    
    ## Returning the function with self starter and names
    returnList <- list(fct = fct, ssfct = ssfct, names = pnames, text = text, noParm = sum(is.na(fixed)))
    
    class(returnList) <- "drcMean"
    invisible(returnList)
  }
