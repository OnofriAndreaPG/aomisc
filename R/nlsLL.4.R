#Log-Logistic Function for bioassay work nlsLL.4
nlsLL.4mean <- function(predictor, b, c, d, ED50) {
                      x <- predictor
                      c+(d-c)/(1+exp(b*(log(x+0.000001)-log(ED50))))
}

#NLSLL.4mean <- deriv(~c+(d-c)/(1+exp(b*(log(predictor+0.000001)-log(ED50)))),c("c","d","b","ED50"),function(predictor,c,d,b,ED50){})

nlsLL.4Init <- function(mCall, LHS, data) {
          xy <- sortedXyData(mCall[["predictor"]], LHS, data)
          x <-  xy[, "x"]; y <- xy[, "y"]
	       c <- min(y)*0.95
          d <- max(y) * 1.05 
              
          ## Linear regression on pseudo y values
          pseudoY <- log((d-y)/(y-c))
          coefs <- coef( lm(pseudoY ~ log(x+0.000001)))
          k <- -coefs[1]; b <- coefs[2]
          ED50 <- exp(k/b)
          value <- c(b,c,d,ED50)
          names(value) <- mCall[c("b", "c", "d", "ED50")]
          value
}

nlsLL.4 <- selfStart(nlsLL.4mean, nlsLL.4Init, parameters=c("b", "c", "d", "ED50"))
