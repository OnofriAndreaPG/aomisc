# Bragg's equation
lorentz.3.fun <- function(X, b, d, e){
  d / ( 1 + b * (X - e)^2)
}

DRC.lorentz.3 <- function(){
  fct <- function(x, parm) {
    lorentz.3.fun(x, parm[,1], parm[,2], parm[,3])
  }
  ssfct <- function(data){
    # Get the data     
    x <- data[, 1]
    y <- data[, 2]
    
    d <- max(y)
    e <- x[which.max(y)]
    
    ## Linear regression on pseudo-y and pseudo-x
    pseudoY <- ( d - y )/ y
    pseudoX <- (x - e)^2
    coefs <- coef( lm(pseudoY ~ pseudoX - 1) )
    b <- coefs[1]
    start <- c(b, d, e)
    return( start )
  }
  names <- c("b", "d", "e")
  text <- "Lorentz equation with three parameters"
    
  ## Returning the function with self starter and names
  returnList <- list(fct = fct, ssfct = ssfct, names = names, text = text)
  class(returnList) <- "drcMean"
  invisible(returnList)
}

lorentz.3.init <- function(mCall, LHS, data) {
    xy <- sortedXyData(mCall[["X"]], LHS, data)
    x <-  xy[, "x"]; y <- xy[, "y"]
    
    d <- max(y)
    e <- x[which.max(y)]

    ## Linear regression on pseudo-y and pseudo-x
    pseudoY <- ( d - y )/ y
    pseudoX <- (x - e)^2
    coefs <- coef( lm(pseudoY ~ pseudoX - 1) )
    b <- coefs[1]
    start <- c(b, d, e)
    names(start) <- mCall[c("b", "d", "e")]
    start
}

NLS.lorentz.3 <- selfStart(lorentz.3.fun, lorentz.3.init, parameters=c("b", "d", "e"))


lorentz.4.fun <- function(X, b, c, d, e){
  c + (d - c) / ( 1 + b * (X - e)^2)
}

DRC.lorentz.4 <- function(){
  fct <- function(x, parm) {
    lorentz.4.fun(x, parm[,1], parm[,2], parm[,3], parm[,4])
  }
  ssfct <- function(data){
    # Get the data     
    x <- data[, 1]
    y <- data[, 2]
    
    d <- max(y)
    c <- min(y) * 0.95
    e <- x[which.max(y)]
    
    
    ## Linear regression on pseudo-y and pseudo-x
    pseudoY <- (d - y)/(y - c)
    pseudoX <- (x - e)^2
    coefs <- coef( lm(pseudoY ~ pseudoX - 1) )
    b <- coefs[1]
    start <- c(b, c, d, e)
    return( start )
  }
  names <- c("b", "c", "d", "e")
  text <- "Lorentz equation with four parameters"
    
  ## Returning the function with self starter and names
  returnList <- list(fct = fct, ssfct = ssfct, names = names, text = text)
  class(returnList) <- "drcMean"
  invisible(returnList)
}

lorentz.4.init <- function(mCall, LHS, data) {
    xy <- sortedXyData(mCall[["X"]], LHS, data)
    x <-  xy[, "x"]; y <- xy[, "y"]
    
    d <- max(y)
    c <- min(y) * 0.95
    e <- x[which.max(y)]
    
    
    ## Linear regression on pseudo-y and pseudo-x
    pseudoY <- (d - y)/(y - c)
    pseudoX <- (x - e)^2
    coefs <- coef( lm(pseudoY ~ pseudoX - 1) )
    b <- coefs[1]
    start <- c(b, c, d, e)
    names(start) <- mCall[c("b", "c", "d", "e")]
    start
}

NLS.lorentz.4 <- selfStart(lorentz.4.fun, lorentz.4.init, parameters=c("b", "c", "d", "e"))
