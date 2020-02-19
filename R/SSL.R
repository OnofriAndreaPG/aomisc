#Log-Logistic Function for bioassay work nlsL.4
L4.fun <- function(predictor, b, c, d, e) {
  x <- predictor
  c + (d - c)/(1 + exp( - b* (x - e)))
}

L4.Init <- function(mCall, LHS, data) {
  xy <- sortedXyData(mCall[["predictor"]], LHS, data)
  x <-  xy[, "x"]; y <- xy[, "y"]
  d <- max(y) * 1.05
  c <- min(y)
  
  ## Linear regression on pseudo y values
  pseudoY <- log((d - y)/(y+0.00001 - c))
  coefs <- coef( lm(pseudoY ~ x))
  k <- coefs[1]; b <- - coefs[2]
  e <- k/b
  value <- c(b,c,d,e)
  names(value) <- mCall[c("b", "c", "d", "e")]
  value
}

NLS.L4 <- selfStart(L4.fun, L4.Init, parameters=c("b", "c", "d", "e"))

#Log-Logistic Function for bioassay work nlsL.3
# Similar to L.3
L3.fun <- function(predictor, b, d, e) {
  x <- predictor
  d/(1 + exp( - b* (x - e)))
}

L3.Init <- function(mCall, LHS, data) {
  xy <- sortedXyData(mCall[["predictor"]], LHS, data)
  x <-  xy[, "x"]; y <- xy[, "y"]
  # y <- beetGrowth$weightFree; x <- beetGrowth$DAE 
  d <- max(y) * 1.05              
  # print(d)
  ## Linear regression on pseudo y values
  pseudoY <- log((d-y)/(y+0.00001))
  coefs <- coef( lm(pseudoY ~ x))
  k <- coefs[1]; b <- -coefs[2]
  e <- k/b
  value <- c(b,d,e)
  names(value) <- mCall[c("b", "d", "e")]
  value
}

NLS.L3 <- selfStart(L3.fun, L3.Init, parameters=c("b", "d", "e"))

#Log-Logistic Function for bioassay work nlsL.2
# Similar to L.3
L2.fun <- function(predictor, b, e) {
  x <- predictor
  1/(1 + exp( - b* (x - e)))
}

L2.Init <- function(mCall, LHS, data) {
  xy <- sortedXyData(mCall[["predictor"]], LHS, data)
  x <-  xy[, "x"]; y <- xy[, "y"]
  d <- 1              
  ## Linear regression on pseudo y values
  pseudoY <- log((d - y)/(y+0.00001))
  coefs <- coef( lm(pseudoY ~ x))
  k <- coefs[1]; b <- - coefs[2]
  e <- k/b
  value <- c(b, e)
  names(value) <- mCall[c("b", "e")]
  value
}

NLS.L2 <- selfStart(L2.fun, L2.Init, parameters=c("b", "e"))

# Logistic curve with two-parameters (DRC) ##########################
"L.2" <-
  function (upper = 1, fixed = c(NA, NA), names = c("b", "e"))
  {
    numParm <- 2
    if (!is.character(names) | !(length(names) == numParm)) {
      stop("Not correct 'names' argument")
    }
    if (!(length(fixed) == numParm)) {
      stop("Not correct length of 'fixed' argument")
    }
    return(logistic(fixed = c(fixed[1], 0, upper, fixed[2], 1),
                    names = c(names[1], "c", "d", names[2], "f"),
                    fctName = as.character(match.call()[[1]]),
                    fctText = "Logistic (ED50 as parameter)"))
  }

