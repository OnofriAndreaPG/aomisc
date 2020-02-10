#Negative exponential - Monomulecolar growth
monoGrowthMean <- function(predictor, a, b, c) {
                      x <- predictor
                      a - (a - b) * exp (- c * x)
}

monoGrowthInit <- function(mCall, LHS, data) {
          xy <- sortedXyData(mCall[["predictor"]], LHS, data)
          x <-  xy[, "x"]; y <- xy[, "y"]
          plateau <- max(y) * 1.05
        
          ## Linear regression on pseudo y values
          pseudoY <- log( 1 - (y / plateau ) )
          coefs <- coef( lm(pseudoY ~ x) )
          temp <- exp(coefs[1])
          b <- plateau * (1 - temp)
          c <- - coefs[2]
          a <- plateau
          value <- c(a, b, c)
          names(value) <- mCall[c("a", "b", "c")]
          value
}

NLSmonoGrowth <- selfStart(monoGrowthMean, monoGrowthInit, parameters=c("a", "b", "c"))

#Logistic growth - 1 ###############################
logiGrowth1Mean <- function(predictor, a, b, c) {
                      x <- predictor
                      a / (1 + exp(- b * x + c))
}

logiGrowth1Init <- function(mCall, LHS, data) {
          xy <- sortedXyData(mCall[["predictor"]], LHS, data)
          x <-  xy[, "x"]; y <- xy[, "y"]
          
          a <- max(y) * 1.05
        
          ## Linear regression on pseudo y values
          pseudoY <- log( (a / y ) - 1 )
          coefs <- coef( lm(pseudoY ~ x) )
          b <- - coefs[2]
          c <- coefs[1]
          value <- c(a, b, c)
          names(value) <- mCall[c("a", "b", "c")]
          value
}

NLSlogiGrowth.1 <- selfStart(logiGrowth1Mean, logiGrowth1Init, parameters=c("a", "b", "c"))


#Logistic Growth 2
logiGrowth2Mean <- function(predictor, a, b, c) {
                      x <- predictor
                      a / (1 + b * exp(- c * x))
}

logiGrowth2Init <- function(mCall, LHS, data) {
          xy <- sortedXyData(mCall[["predictor"]], LHS, data)
          x <-  xy[, "x"]; y <- xy[, "y"]
          
          a <- max(y) * 1.05
        
          ## Linear regression on pseudo y values
          pseudoY <- log( (a / (y + 0.0001) ) - 1 )
          coefs <- coef( lm(pseudoY ~ x) )
          c <- - coefs[2]
          k <- coefs[1]
          b <- exp(k)
          value <- c(a, b, c)
          names(value) <- mCall[c("a", "b", "c")]
          value
}

NLSlogiGrowth.2 <- selfStart(logiGrowth2Mean, logiGrowth2Init, parameters=c("a", "b", "c"))

#Logistic Growth 3
logiGrowth3Mean <- function(predictor, init, m, plateau) {
                      x <- predictor
                      init * plateau / (init + (plateau - init) * exp( - m * x))
}

logiGrowth3Init <- function(mCall, LHS, data) {
			xy <- sortedXyData(mCall[["predictor"]], LHS, data)
         x <-  xy[, "x"]; y <- xy[, "y"]
          
        	plateau <- max(y) * 1.05
        
        	## Linear regression on pseudo y values
        	pseudoY <- log( (plateau / (y + 0.0001) ) - 1 )
        	coefs <- coef( lm(pseudoY ~ x) )
        	b <- exp(coefs[1])
        	init <- plateau / (1 + b)
        	m <- - coefs[2]
         value <- c(init, m, plateau)
         names(value) <- mCall[c("init", "m", "plateau")]
         value
}

NLSlogiGrowth.3 <- selfStart(logiGrowth3Mean, logiGrowth3Init, parameters=c("init", "m", "plateau"))

#Logistic Growth 4
logiGrowth4Mean <- function(predictor, t50, m, plateau) {
                      x <- predictor
                      plateau / (1 + exp(- m * (x - t50)))}

logiGrowth4Init <- function(mCall, LHS, data) {
			xy <- sortedXyData(mCall[["predictor"]], LHS, data)
         x <-  xy[, "x"]; y <- xy[, "y"]
          
        	plateau <- max(y) * 1.05
        
         ## Linear regression on pseudo y values
         pseudoY <- log( (plateau / (y + 0.0001) ) - 1 )
         coefs <- coef( lm(pseudoY ~ x) )
         m <- - coefs[2]
         t50 <- coefs[1] / m
         value <- c(t50, m, plateau)
         names(value) <- mCall[c("t50", "m", "plateau")]
         value
}

NLSlogiGrowth.4 <- selfStart(logiGrowth4Mean, logiGrowth4Init, parameters=c("t50", "m", "plateau"))

#Logistic function 5
logistic5Mean <- function(predictor, a, b, c, d) {
                      x <- predictor
                      c + (d-c) / (1 + b * exp(a * x))
}

logistic5Init <- function(mCall, LHS, data) {
          xy <- sortedXyData(mCall[["predictor"]], LHS, data)
          x <-  xy[, "x"]; y <- xy[, "y"]
          
          d <- max(y) * 1.05
          c <- min(y) -0.0001
          ## Linear regression on pseudo y values
          pseudoY <- log((d - y) / (y - c ))
          coefs <- coef( lm(pseudoY ~ x) )
          a <- coefs[2]
          k <- coefs[1]
          b <- exp(k)
          value <- c(a, b, c, d)
          names(value) <- mCall[c("a", "b", "c", "d")]
          value
}

NLSlogistic.1 <- selfStart(logistic5Mean, logistic5Init, parameters=c("a", "b", "c", "d"))

#Logistic function 6
logistic6Mean <- function(predictor, a, b) {
                      x <- predictor
                      1 / (1 + b * exp(- a * x))
}

logistic6Init <- function(mCall, LHS, data) {
          xy <- sortedXyData(mCall[["predictor"]], LHS, data)
          x <-  xy[, "x"]; y <- xy[, "y"]
          
          ## Linear regression on pseudo y values
          pseudoY <- log((1 - y) / y)
          coefs <- coef( lm(pseudoY ~ x) )
          a <- - coefs[2]
          k <- coefs[1]
          b <- exp(k)
          value <- c(a, b)
          names(value) <- mCall[c("a", "b")]
          value
}

NLSlogistic.2 <- selfStart(logistic6Mean, logistic6Init, parameters=c("a", "b"))

#LOG_LOGISTIC FUNCTIONS##########################################################

#Hill function
hillCurveMean <- function(predictor, a, b, c) {
                      (a * predictor^c)/(b + predictor^c)
}

hillCurveInit <- function(mCall, LHS, data) {
          xy <- sortedXyData(mCall[["predictor"]], LHS, data)
	  		 x <-  xy[, "x"]; y <- xy[, "y"]
          a <- max(y) * 1.05
          pseudoY <-  log(( a - y )/ y)
          pseudoX <- log(x)
          lmFit <- lm(pseudoY ~ pseudoX )
          coefs <- coef(lmFit)
          b <- exp(coefs[1])
          c <- - coefs[2]
          value <- c(a, b, c)
          names(value) <- mCall[c("a", "b", "c")]
          value
}

NLShillCurve <- selfStart(hillCurveMean, hillCurveInit, parameters=c("a", "b", "c"))

#GOMPERTZ MODELS################################################################

#Gompertz growth - 1
gompGrowth1Mean <- function(predictor, a, m, c) {
                      x <- predictor
                      a * exp( - (m/c) * exp (-c * x))
}

gompGrowth1Init <- function(mCall, LHS, data) {
          xy <- sortedXyData(mCall[["predictor"]], LHS, data)
          x <-  xy[, "x"]; y <- xy[, "y"]
          
          plateau <- max(y) * 1.05
        
        	 ## Linear regression on pseudo y values
          pseudoY <- log( - log( y / plateau ) )
          coefs <- coef( lm(pseudoY ~ x) )
          k <- coefs[1]; c <- - coefs[2]
          b <- exp(k) 
          m <- b * c
          a <- plateau
          value <- c(a, m, c)
          names(value) <- mCall[c("a", "m", "c")]
          value
}

NLSgompGrowth.1 <- selfStart(gompGrowth1Mean, gompGrowth1Init, parameters=c("a", "m", "c"))

#Gompertz growth 2
gompGrowth2Mean <- function(predictor, a, b, c) {
                      x <- predictor
                      a * exp( - exp (b - c*x))
}

gompGrowth2Init <- function(mCall, LHS, data) {
          xy <- sortedXyData(mCall[["predictor"]], LHS, data)
          x <-  xy[, "x"]; y <- xy[, "y"]
          
          a <- max(y) * 1.05
        
          ## Linear regression on pseudo y values
          pseudoY <- log( - log( y / a ) )
          coefs <- coef( lm(pseudoY ~ x) )

          k <- coefs[1]
          c <- - coefs[2]
          b <- k
          value <- c(a, b, c)
          names(value) <- mCall[c("a", "b", "c")]
          value
}

NLSgompGrowth.2 <- selfStart(gompGrowth2Mean, gompGrowth2Init, parameters=c("a", "b", "c"))

#Gompertz growth - 3
gompGrowth3Mean <- function(predictor, a, b, c) {
                      x <- predictor
                      a * exp( -b * exp (-c*x))
}

gompGrowth3Init <- function(mCall, LHS, data) {
          xy <- sortedXyData(mCall[["predictor"]], LHS, data)
          x <-  xy[, "x"]; y <- xy[, "y"]
          
          a <- max(y) * 1.05
        
          ## Linear regression on pseudo y values
          pseudoY <- log( - log( y / a ) )
          coefs <- coef( lm(pseudoY ~ x) )

          k <- coefs[1]
          c <- - coefs[2]
          b <- exp(k)
          value <- c(a, b, c)
          names(value) <- mCall[c("a", "b", "c")]
          value
}

NLSgompGrowth.3 <- selfStart(gompGrowth3Mean, gompGrowth3Init, parameters=c("a", "b", "c"))


#Extreme value
extremeValueMean <- function(predictor, a, b, c) {
                      x <- predictor
                      a * (1 - exp( - exp (b - c*x)))
}

extremeValueInit <- function(mCall, LHS, data) {
          xy <- sortedXyData(mCall[["predictor"]], LHS, data)
          x <-  xy[, "x"]; y <- xy[, "y"]
          
          a <- max(y) * 1.05
        
          ## Linear regression on pseudo y values
          pseudoY <- log( - log( (a - y ) / a ) )
          coefs <- coef( lm(pseudoY ~ x) )

          k <- coefs[1]
          c <- - coefs[2]
          b <- k
          value <- c(a, b, c)
          names(value) <- mCall[c("a", "b", "c")]
          value
}

NLSextremeValue <- selfStart(extremeValueMean, extremeValueInit, parameters=c("a", "b", "c"))

#WEIBULL TYPE MODELS
#Weibull-1

#Modified Mitscherlich equation for A vs PFD relationships
AvsPFDMean <- function(predictor, Rd, Amax, Qapp) {
                      x <- predictor
                      Rd+Amax*(1-exp((-Qapp/Amax)*x))
}

AvsPFDInit <- function(mCall, LHS, data) {
          xy <- sortedXyData(mCall[["predictor"]], LHS, data)
          x <-  xy[, "x"]; y <- xy[, "y"]
	  interc <- min(y)
          plateau <- max(y) * 1.05 - interc        
          ## Linear regression on pseudo y values
          pseudoY <- log( 1 - (((y - interc) / plateau ) ))
          coefs <- coef( lm(pseudoY ~ x - 1) )
          Amax <- plateau
          Rd <- interc
          b <-  coefs[1]
          Qapp <- -b*plateau
          value <- c(Rd,Amax,Qapp)
          names(value) <- mCall[c("Rd", "Amax", "Qapp")]
          value
}

NLSAvsPFD <- selfStart(AvsPFDMean, AvsPFDInit, parameters=c("Rd", "Amax", "Qapp"))

#Inverse polynomial
polyInv.3mean <- function(predictor, a, b, c) {
                      x <- predictor
                      1/(a + b*x + c*x^2)
}

polyInv.3Init <- function(mCall, LHS, data) {
          xy <- sortedXyData(mCall[["predictor"]], LHS, data)
          x <-  xy[, "x"]; y <- xy[, "y"]
          ## Linear regression on pseudo y values
          #pseudoY <- 1/y
          coefs <- coef(glm(y ~ x + I(x^2), family=gaussian(link="inverse")))
          a <- coefs[1]; b <- coefs[2]; c <- coefs[3]
          
          value <- c(a,b,c)
          names(value) <- mCall[c("a", "b", "c")]
          value
}

NLSpolyInv.3 <- selfStart(polyInv.3mean, polyInv.3Init, parameters=c("a", "b", "c"))

#Inverse polynomial 2
polyInv.4mean <- function(predictor, a, b, c) {
                      x <- predictor
                      x/(a + b*x + c*x^2)
}

polyInv.4Init <- function(mCall, LHS, data) {
          xy <- sortedXyData(mCall[["predictor"]], LHS, data)
          x <-  xy[, "x"]; y <- xy[, "y"]
          ## Linear regression on pseudo y values
          pseudoY <- y/x
          coefs <- coef(glm(pseudoY ~ x + I(x^2), family=gaussian(link="inverse")))
          a <- coefs[1]; b <- coefs[2]; c <- coefs[3]
          print(a);print(b);print(c)
          value <- c(a,b,c)
          names(value) <- mCall[c("a", "b", "c")]
          value
}

NLSpolyInv.4 <- selfStart(polyInv.4mean, polyInv.4Init, parameters=c("a", "b", "c"))

