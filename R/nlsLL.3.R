#Log-Logistic Function for bioassay work nlsLL.3
nlsLL.3mean <- function(predictor, b, d, ED50) {
                      x <- predictor
                      d/(1+exp(b*(log(x+0.000001)-log(ED50))))
}

nlsLL.3Init <- function(mCall, LHS, data) {
          xy <- sortedXyData(mCall[["predictor"]], LHS, data)
          x <-  xy[, "x"]; y <- xy[, "y"]
          d <- max(y) * 1.05              
          ## Linear regression on pseudo y values
          pseudoY <- log((d-y)/(y+0.00001))
          coefs <- coef( lm(pseudoY ~ log(x+0.000001)))
          k <- -coefs[1]; b <- coefs[2]
          ED50 <- exp(k/b)
          value <- c(b,d,ED50)
          names(value) <- mCall[c("b", "d", "ED50")]
          value
}

nlsLL.3 <- selfStart(nlsLL.3mean, nlsLL.3Init, parameters=c("b", "d", "ED50"))
