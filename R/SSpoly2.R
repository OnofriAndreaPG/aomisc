#Polynomial regression (2nd order) ##############################################
poly2.fun <- function(predictor, a, b, c) {
  a + b * predictor + c * (predictor^2)
}

poly2.Init <- function(mCall, LHS, data) {
  xy <- sortedXyData(mCall[["predictor"]], LHS, data)
  lmFit <- lm((xy[, "y"]) ~ xy[, "x"] + I( xy[, "x"]^2))
  coefs <- coef(lmFit)
  a <- coefs[1]
  b <- coefs[2]
  c <- coefs[3]
  value <- c(a, b, c)
  names(value) <- mCall[c("a", "b", "c")]
  value
}
NLS.poly2 <- selfStart(poly2.fun, poly2.Init, parameters=c("a", "b", "c"))

"DRC.poly2" <- function(fixed = c(NA, NA, NA), names = c("a", "b", "c"))
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
    a + b * x + c * (x^2)
  }
  
  ## Defining self starter function
  ssfct <- function(dataf)
  {
    x <- dataf[, 1]
    y <- dataf[, 2]
    
    #regression on pseudo y values
    pseudoY <- y
    pseudoX <- x
    coefs <- coef( lm(pseudoY ~ pseudoX + I(pseudoX^2)) )
    a <- coefs[1]
    b <- coefs[2]
    c <- coefs[3]
    
    return(c(a, b, c)[notFixed])
  }
  
  ## Defining names
  pnames <- names[notFixed]
  
  ## Defining derivatives
  deriv1 <- function(x, parm){
    d1 <- rep(1, length(x) )
    d2 <- x
    d3 <- x^2
    cbind(d1, d2, d3)
  }
  ## Defining the ED function
  
  ## Defining the inverse function
  
  ## Defining descriptive text
  text <- "Second Order Polynomial"
  
  ## Returning the function with self starter and names
  returnList <- list(fct = fct, ssfct = ssfct, names = pnames, 
                     text = text, noParm = sum(is.na(fixed)), deriv1 = deriv1)
  
  class(returnList) <- "drcMean"
  invisible(returnList)
}
