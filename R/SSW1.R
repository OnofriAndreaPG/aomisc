# Weibul type 1 Function for bioassay work
# Edited on 07/02/2020 
W1.4.fun <- function(predictor, b, c, d, e) {
  x <- predictor
  c + (d - c) * exp ( - exp ( - b*(log(x + 0.000001) - log(e))))
}

W1.4.init <- function(mCall, LHS, data) {
  xy <- sortedXyData(mCall[["predictor"]], LHS, data)
  x <-  xy[, "x"]; y <- xy[, "y"]
  
  d <- max(y) * 1.05
  c <- min(y) * 0.95
  
  ## Linear regression on pseudo y values
  pseudoY <- log( - log( (y - c) / (d - c) ) )
  coefs <- coef( lm(pseudoY ~ log(x+0.000001)) )
  
  b <- - coefs[2]
  e <- exp(coefs[1]/b)
  value <- c(b, c, d, e)
  names(value) <- mCall[c("b", "c", "d", "e")]
  value
}

NLS.W1.4 <- selfStart(W1.4.fun, W1.4.init, parameters=c("b", "c", "d", "e"))


# Weibul type 2 Function for bioassay work nlsW1.3
# Edited on 07/02/2020 
W1.3.fun <- function(predictor, b, d, e) {
                      x <- predictor
                      d * exp ( - exp ( - b*(log(x+0.000001)-log(e))))
}

W1.3.init <- function(mCall, LHS, data) {
          xy <- sortedXyData(mCall[["predictor"]], LHS, data)
          x <-  xy[, "x"]; y <- xy[, "y"]
          
          d <- max(y) * 1.05
        
          ## Linear regression on pseudo y values
          pseudoY <- log( - log( y / d ) )
          coefs <- coef( lm(pseudoY ~ log(x+0.000001)) )

          b <- - coefs[2]
          e <- exp(coefs[1]/b)
          value <- c(b, d, e)
          names(value) <- mCall[c("b", "d", "e")]
          value
}

NLS.W1.3 <- selfStart(W1.3.fun, W1.3.init, parameters=c("b", "d", "e"))

# Weibul type 2 Function for bioassay work nlsW1.3
# Edited on 07/02/2020 
W1.2.fun <- function(predictor, b, e) {
  x <- predictor
  exp ( - exp ( - b*(log(x+0.000001)-log(e))))
}

W1.2.init <- function(mCall, LHS, data) {
  xy <- sortedXyData(mCall[["predictor"]], LHS, data)
  x <-  xy[, "x"]; y <- xy[, "y"]
  
  d <- 1
  
  ## Linear regression on pseudo y values
  pseudoY <- log( - log( y / d ) )
  coefs <- coef( lm(pseudoY ~ log(x+0.000001)) )
  
  b <- - coefs[2]
  e <- exp(coefs[1]/b)
  value <- c(b, e)
  names(value) <- mCall[c("b", "e")]
  value
}

NLS.W1.2 <- selfStart(W1.2.fun, W1.2.init, parameters=c("b", "e"))

