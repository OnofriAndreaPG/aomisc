# Weibul type 1 Function for bioassay work nlsW2.4
# Edited on 07/02/2020 
W2.4.fun <- function(predictor, b, c, d, e) {
  x <- predictor
  c + (d - c) * (1 - exp( - exp (b * (log(x + 0.0000001) - log(e)))))
}

W2.4.init <- function(mCall, LHS, data) {
  xy <- sortedXyData(mCall[["predictor"]], LHS, data)
  x <-  xy[, "x"]; y <- xy[, "y"]
  
  x <- brassica$Dose
  y <- brassica$FW
  d <- max(y) * 1.05
  c <- min(y) * 0.90
  
  ## Linear regression on pseudo y values
  pseudoY <- log( - log( (d - y ) / (d - c) ) )
  coefs <- coef( lm(pseudoY ~ log(x+0.0000001)) )
  
  b <- coefs[2]
  e <- exp( - coefs[1]/b)

  value <- c(b, c, d, e)
  names(value) <- mCall[c("b", "c", "d", "e")]
  value
}

NLS.W2.4 <- selfStart(W2.4.fun, W2.4.init, parameters=c("b", "c", "d", "e"))

# Weibul type 1 Function for bioassay work nlsW2.3
# Edited on 07/02/2020 
W2.3.fun <- function(predictor, b, d, e) {
                      x <- predictor
                      d * (1 - exp( - exp (b * (log(x + 0.0000001) - log(e)))))
}

W2.3.init <- function(mCall, LHS, data) {
          xy <- sortedXyData(mCall[["predictor"]], LHS, data)
          x <-  xy[, "x"]; y <- xy[, "y"]
          
          d <- max(y) * 1.05
        
          ## Linear regression on pseudo y values
          pseudoY <- log( - log( (d - y ) / d ) )
          coefs <- coef( lm(pseudoY ~ log(x+0.0000001)) )

          b <- coefs[2]
          e <- exp( - coefs[1]/b)
          value <- c(b, d, e)
          names(value) <- mCall[c("b", "d", "e")]
          value
}

NLS.W2.3 <- selfStart(W2.3.fun, W2.3.init, parameters=c("b", "d", "e"))

# Weibul type 1 Function for bioassay work nlsW2.3
# Edited on 07/02/2020 
W2.2.fun <- function(predictor, b, e) {
  x <- predictor
  1 - exp( - exp (b * (log(x + 0.0000001) - log(e))))
}

W2.2.init <- function(mCall, LHS, data) {
  xy <- sortedXyData(mCall[["predictor"]], LHS, data)
  x <-  xy[, "x"]; y <- xy[, "y"]
  
  d <- 1
  
  ## Linear regression on pseudo y values
  pseudoY <- log( - log( (d - y ) / d ) )
  coefs <- coef( lm(pseudoY ~ log(x+0.0000001)) )
  
  b <- coefs[2]
  e <- exp( - coefs[1]/b)

  value <- c(b, e)
  names(value) <- mCall[c("b", "e")]
  value
}

NLS.W2.2 <- selfStart(W2.2.fun, W2.2.init, parameters=c("b", "e"))

