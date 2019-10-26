#Log-Logistic Function for bioassay work nlsL.3
# Similar to L.3
nlsL.3mean <- function(predictor, b, d, e) {
  x <- predictor
  d/(1+exp(b*(x - e)))
}

nlsL.3Init <- function(mCall, LHS, data) {
  xy <- sortedXyData(mCall[["predictor"]], LHS, data)
  x <-  xy[, "x"]; y <- xy[, "y"]
  d <- max(y) * 1.05              
  ## Linear regression on pseudo y values
  pseudoY <- log((d-y)/(y+0.00001))
  coefs <- coef( lm(pseudoY ~ x))
  k <- -coefs[1]; b <- coefs[2]
  e <- exp(k/b)
  value <- c(b,d,e)
  names(value) <- mCall[c("b", "d", "e")]
  value
}

nlsL.3 <- selfStart(nlsL.3mean, nlsL.3Init, parameters=c("b", "d", "e"))
