# Gompertz equation for bioassay work 
G4.fun <- function(predictor, b, c, d, e) {
  x <- predictor
  c + (d - c) * (exp (- exp ( - b * ( x - e))))
}

G4.Init <- function(mCall, LHS, data) {
  xy <- sortedXyData(mCall[["predictor"]], LHS, data)
  x <-  xy[, "x"]; y <- xy[, "y"]
  d <- max(y) * 1.05
  c <- min(y) * 0.95
  
  ## Linear regression on pseudo y values
  pseudoY <- log(-log((y - c)/(d - c)))
  coefs <- coef( lm(pseudoY ~ x))
  k <- coefs[1]; b <- - coefs[2]
  e <- k/b
  value <- c(b,c,d,e)
  names(value) <- mCall[c("b", "c", "d", "e")]
  value
}

NLS.G4 <- selfStart(G4.fun, G4.Init, parameters=c("b", "c", "d", "e"))

G3.fun <- function(predictor, b, d, e) {
  x <- predictor
  d * (exp (- exp ( - b * ( x - e))))
}

G3.Init <- function(mCall, LHS, data) {
  xy <- sortedXyData(mCall[["predictor"]], LHS, data)
  x <-  xy[, "x"]; y <- xy[, "y"]
  # y <- beetGrowth$weightFree; x <- beetGrowth$DAE 
  d <- max(y) * 1.05              
  # print(d)
  ## Linear regression on pseudo y values
  pseudoY <- log(-log(y/d))
  coefs <- coef( lm(pseudoY ~ x))
  k <- coefs[1]; b <- -coefs[2]
  e <- k/b
  value <- c(b,d,e)
  names(value) <- mCall[c("b", "d", "e")]
  value
}

NLS.G3 <- selfStart(G3.fun, G3.Init, parameters=c("b", "d", "e"))

G2.fun <- function(predictor, b, e) {
  x <- predictor
  exp (- exp ( - b * ( x - e)))
}

G2.Init <- function(mCall, LHS, data) {
  xy <- sortedXyData(mCall[["predictor"]], LHS, data)
  x <-  xy[, "x"]; y <- xy[, "y"]
                
  ## Linear regression on pseudo y values
  pseudoY <- log(-log(y) )
  coefs <- coef( lm(pseudoY ~ x))
  k <- coefs[1]; b <- - coefs[2]
  e <- k/b
  value <- c(b, e)
  names(value) <- mCall[c("b", "e")]
  value
}

NLS.G2 <- selfStart(G2.fun, G2.Init, parameters=c("b", "e"))
