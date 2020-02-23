beta.fun <- function(X, b, d, Xb, Xo, Xc){
  .expr1 <-  (X - Xb)/(Xo - Xb)
  .expr2 <- (Xc - X)/(Xc - Xo)
  .expr3 <- (Xc - Xo)/(Xo - Xb)
  #ifelse(temp > tb & temp < tc, (.expr2*.expr1^(1/.expr3))^a, 0)
  ifelse(X > Xb & X < Xc, d * (.expr1*.expr2^.expr3)^b, 0)
}

DRC.beta <- function(){

  fct <- function(x, parm) {
      # function code here
    beta.fun(x, parm[,1], parm[,2], parm[,3], parm[,4], parm[,5])
  }
  ssfct <- function(data){
     # Self-starting code here
    x <- data[, 1]
    y <- data[, 2]

    d <- max(y)
    Xo <- x[which.max(y)]
    firstidx <- min( which(y !=0) )
    Xb <- ifelse(firstidx == 1,  x[1], (x[firstidx] + x[(firstidx - 1)])/2)
    secidx <- max( which(y !=0) )
    Xc <- ifelse(secidx == length(y),  x[length(x)], (x[secidx] + x[(secidx + 1)])/2)
    c(1, d, Xb, Xo, Xc)
  
    }
  names <- c("b", "d", "Xb", "Xo", "Xc")
  text <- "Beta function"
    
  ## Returning the function with self starter and names
  returnList <- list(fct = fct, ssfct = ssfct, names = names, text = text)
  class(returnList) <- "drcMean"
  invisible(returnList)
}

beta.init <- function(mCall, LHS, data) {
    xy <- sortedXyData(mCall[["X"]], LHS, data)
    x <-  xy[, "x"]; y <- xy[, "y"]
    #Self starting code ##############
    d <- max(y)
    Xo <- x[which.max(y)]
    firstidx <- min( which(y !=0) )
    Xb <- ifelse(firstidx == 1,  x[1], (x[firstidx] + x[(firstidx - 1)])/2)
    secidx <- max( which(y !=0) )
    Xc <- ifelse(secidx == length(y),  x[length(x)], (x[secidx] + x[(secidx + 1)])/2)

    start <- c(1, d, Xb, Xo, Xc)
    names(start) <- mCall[c("b", "d", "Xb", "Xo", "Xc")]
    start
}

NLS.beta <- selfStart(beta.fun, beta.init, parameters=c("b", "d", "Xb", "Xo", "Xc"))

