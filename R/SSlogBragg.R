# Bragg's equation
logBragg.3.fun <- function(X, b, d, e){
  d * exp(- b * (log(X + 0.000001) - e)^2)
}

# Da fare
DRC.logBragg.3 <- function(){
  fct <- function(x, parm) {
    logBragg.3.fun(x, parm[,1], parm[,2], parm[,3])
  }
  ssfct <- function(data){
    # Get the data     
    x <- log(data[, 1] + 0.000001)
    y <- data[, 2]
    
    d <- max(y)
    e <- x[which.max(y)]
    
    ## Linear regression on pseudo-y and pseudo-x
    pseudoY <- log( (y + 0.0001) / d )
    pseudoX <- (x - e)^2
    coefs <- coef( lm(pseudoY ~ pseudoX - 1) )
    b <- - coefs[1]
    start <- c(b, d, e)
    return( start )
  }
  names <- c("b", "d", "e")
  text <- "log-Bragg equation with three parameters"
    
  ## Returning the function with self starter and names
  returnList <- list(fct = fct, ssfct = ssfct, names = names, text = text)
  class(returnList) <- "drcMean"
  invisible(returnList)
}

bragg.3.init <- function(mCall, LHS, data, ...) {
    xy <- sortedXyData(mCall[["X"]], LHS, data)
    x <-  xy[, "x"]; y <- xy[, "y"]
    
    d <- max(y)
    e <- x[which.max(y)]

    ## Linear regression on pseudo-y and pseudo-x
    pseudoY <- log( (y + 0.0001) / d )
    pseudoX <- (x - e)^2
    coefs <- coef( lm(pseudoY ~ pseudoX - 1) )
    b <- - coefs[1]
    start <- c(b, d, e)
    names(start) <- mCall[c("b", "d", "e")]
    start
}

NLS.bragg.3 <- selfStart(bragg.3.fun, bragg.3.init, parameters=c("b", "d", "e"))


bragg.4.fun <- function(X, b, c, d, e){
  c + (d - c) * exp(- b * (X - e)^2)
}

DRC.bragg.4 <- function(){
  fct <- function(x, parm) {
    bragg.4.fun(x, parm[,1], parm[,2], parm[,3], parm[,4])
  }
  ssfct <- function(data){
    # Get the data     
    x <- data[, 1]
    y <- data[, 2]
    
    d <- max(y)
    c <- min(y) * 0.95
    e <- x[which.max(y)]
    
    
    ## Linear regression on pseudo-y and pseudo-x
    pseudoY <- log( ((y + 0.0001) - c) / d )
    pseudoX <- (x - e)^2
    coefs <- coef( lm(pseudoY ~ pseudoX - 1) )
    b <- - coefs[1]
    start <- c(b, c, d, e)
    return( start )
  }
  names <- c("b", "c", "d", "e")
  text <- "Bragg equation with four parameters"
    
  ## Returning the function with self starter and names
  returnList <- list(fct = fct, ssfct = ssfct, names = names, text = text)
  class(returnList) <- "drcMean"
  invisible(returnList)
}

bragg.4.init <- function(mCall, LHS, data, ...) {
    xy <- sortedXyData(mCall[["X"]], LHS, data)
    x <-  xy[, "x"]; y <- xy[, "y"]
    
    d <- max(y)
    c <- min(y) * 0.95
    e <- x[which.max(y)]
    
    
    ## Linear regression on pseudo-y and pseudo-x
    pseudoY <- log( ((y + 0.0001) - c) / d )
    pseudoX <- (x - e)^2
    coefs <- coef( lm(pseudoY ~ pseudoX - 1) )
    b <- - coefs[1]
    start <- c(b, c, d, e)
    names(start) <- mCall[c("b", "c", "d", "e")]
    start
}

NLS.bragg.4 <- selfStart(bragg.4.fun, bragg.4.init, parameters=c("b", "c", "d", "e"))
