#Log-Logistic Function for bioassay work nlsLL.4
LL.4.fun <- function(predictor, b, c, d, e) {
  x <- predictor
  c+(d-c)/(1+exp(- b*(log(x+0.000001)-log(e))))
}

#NLSLL.4mean <- deriv(~c+(d-c)/(1+exp(b*(log(predictor+0.000001)-log(ED50)))),c("c","d","b","ED50"),function(predictor,c,d,b,ED50){})

LL.4.Init <- function(mCall, LHS, data) {
  xy <- sortedXyData(mCall[["predictor"]], LHS, data)
  x <-  xy[, "x"]; y <- xy[, "y"]
  c <- min(y)*0.95
  d <- max(y) * 1.05 
  
  ## Linear regression on pseudo y values
  pseudoY <- log((d-y)/(y-c))
  coefs <- coef( lm(pseudoY ~ log(x+0.000001)))
  k <- -coefs[1]; b <- -coefs[2]
  e <- exp(k/b)
  value <- c(b,c,d,e)
  names(value) <- mCall[c("b", "c", "d", "e")]
  value
}

nlsLL.4 <- selfStart(LL.4.fun, LL.4.Init, parameters=c("b", "c", "d", "e"))

# Log-Logistic Function for bioassay work nlsLL.3
# Edited on 07/02/2020 
LL.3.fun <- function(predictor, b, d, e) {
                      x <- predictor
                      d/(1+exp(-b*(log(x+0.000001)-log(e))))
}

LL.3.init <- function(mCall, LHS, data) {
          xy <- sortedXyData(mCall[["predictor"]], LHS, data)
          x <-  xy[, "x"]; y <- xy[, "y"]
          d <- max(y) * 1.05              
          ## Linear regression on pseudo y values
          pseudoY <- log((d-y)/(y+0.00001))
          coefs <- coef( lm(pseudoY ~ log(x+0.000001)))
          k <- -coefs[1]; b <- -coefs[2]
          e <- exp(k/b)
          value <- c(b,d,e)
          names(value) <- mCall[c("b", "d", "e")]
          value
}

nlsLL.3 <- selfStart(LL.3.fun, LL.3.init, parameters=c("b", "d", "e"))

# Log-Logistic Function for bioassay work nlsLL.2
# Edited on 07/02/2020 
LL.2.fun <- function(predictor, b, e) {
  x <- predictor
  1/(1+exp(-b*(log(x+0.000001)-log(e))))
}

LL.2.init <- function(mCall, LHS, data) {
  xy <- sortedXyData(mCall[["predictor"]], LHS, data)
  x <-  xy[, "x"]; y <- xy[, "y"]
  d <- 1              
  ## Linear regression on pseudo y values
  pseudoY <- log((d-y)/(y+0.00001))
  coefs <- coef( lm(pseudoY ~ log(x+0.000001)))
  k <- -coefs[1]; b <- -coefs[2]
  e <- exp(k/b)
  value <- c(b,e)
  names(value) <- mCall[c("b", "e")]
  value
}

nlsLL.2 <- selfStart(LL.2.fun, LL.2.init, parameters=c("b", "e"))
