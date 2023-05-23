#Linear Model ##############################################
linear.fun <- function(predictor, a, b) {
  a + b * predictor
}

linear.Init <- function(mCall, LHS, data, ...) {
  xy <- sortedXyData(mCall[["predictor"]], LHS, data)
  lmFit <- lm( (xy[, "y"]) ~ xy[, "x"] )
  coefs <- coef(lmFit)
  a <- coefs[1]
  b <- coefs[2]
  value <- c(a, b)
  names(value) <- mCall[c("a", "b")]
  value
}

NLS.linear <- selfStart(linear.fun, linear.Init, parameters=c("a", "b"))

"DRC.linear" <- function(fixed = c(NA, NA), 
                         names = c("a", "b"))
{
  ## Checking arguments
  numParm <- 2
  if (!is.character(names) | !(length(names) == numParm))
  {stop("Not correct 'names' argument")}
  if (!(length(fixed) == numParm))
  {stop("Not correct 'fixed' argument")}
  
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
    linear.fun(x, a, b)
  }
  
  ## Defining self starter function
  ssfct <- function(dataf)
  {
    x <- dataf[, 1]
    y <- dataf[, 2]
    
    #regression on pseudo y values
    pseudoY <- y
    pseudoX <- x
    coefs <- coef( lm(pseudoY ~ pseudoX) )
    a <- coefs[1]
    
    b <- coefs[2]
    
    return(c(a, b)[notFixed])
  }
  
  ## Defining names
  pnames <- names[notFixed]
  
  ## Defining derivatives
  deriv1 <- function(x, parm){
    d1 <- rep(1, length(x) )
    d2 <- x
    cbind(d1, d2)
  }
  ## Defining the ED function
  
  ## Defining the inverse function
  
  ## Defining descriptive text
  text <- "Straight line"
  
  ## Returning the function with self starter and names
  returnList <- list(fct = fct, ssfct = ssfct, deriv1 = deriv1, names = pnames, text = text, noParm = sum(is.na(fixed)))
  
  class(returnList) <- "drcMean"
  invisible(returnList)
}

#linear model with no intercept#########################################
linearOrigin.fun <- function(predictor, b) {
  b * predictor
}

linearOrigin.Init <- function(mCall, LHS, data, ...) {
  xy <- sortedXyData(mCall[["predictor"]], LHS, data)
  lmFit <- lm((xy[, "y"]) ~ xy[, "x"]-1)
  coefs <- coef(lmFit)
  b <- coefs[1]
  value <- c(b)
  names(value) <- mCall[c("b")]
  value
}

NLS.linearOrigin <- selfStart(linearOrigin.fun, linearOrigin.Init, parameters=c("b"))

