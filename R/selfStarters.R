# A collection of self-starting functions for R ########################
# Edited: 9/1/19
# LIST of available functions
# 1 - linear
# 2 - linear by origin
# 3 - second order polynomial
# 4 - exponential growth
# 5 - Exponential decay
# 6 - Asymptotic regression
# 7 - Negative exponential
# 8 - Negative exponential distribution
# 9 - Asymptotic regression
# 9 - Power curve
# 10 - Logarithmic regression
# 11 - Logarithmic regression with no intercept
# 12 - Yield loss function 
# 13 - Yield weed density curve
# 14 - Rational function


#Linear Model ##############################################
linear.fun <- function(predictor, a, b) {
                      a + b * predictor
}

linear.Init <- function(mCall, LHS, data) {
          xy <- sortedXyData(mCall[["predictor"]], LHS, data)
          lmFit <- lm( (xy[, "y"]) ~ xy[, "x"] )
          coefs <- coef(lmFit)
          a <- coefs[1]
          b <- coefs[2]
          value <- c(a, b)
          names(value) <- mCall[c("a", "b")]
          value
}

NLS.linear <- selfStart(linear.fun, linear.Init, parameters=c("a", "b"))

"DRC.linear" <- function(fixed = c(NA, NA), names = c("a", "b"))
{
    ## Checking arguments
    numParm <- 2
    if (!is.character(names) | !(length(names) == numParm))
      {stop("Not correct 'names' argument")}
    if (!(length(fixed) == numParm))
      {stop("Not correct 'fixed' argument")}

    ## Fixing parameters (using argument 'fixed')
    notFixed <- is.na(fixed)
    parmVec <- rep(0, numParm)
    parmVec[!notFixed] <- fixed[!notFixed]

    ## Defining the non-linear function
    fct <- function(x, parm)
    {
        parmMat <- matrix(parmVec, nrow(parm), numParm, byrow = TRUE)
        parmMat[, notFixed] <- parm

        a <- parmMat[, 1]; b <- parmMat[, 2]
        linear.fun(x, a, b)
    }

    ## Defining self starter function
    ssfct <- function(dataf)
    {
        x <- dataf[, 1]
        y <- dataf[, 2]

        #regression on pseudo y values
        pseudoY <- y
        pseudoX <- x
        coefs <- coef( lm(pseudoY ~ pseudoX) )
        a <- coefs[1]

        b <- coefs[2]

        return(c(a, b)[notFixed])
    }

    ## Defining names
    pnames <- names[notFixed]

    ## Defining derivatives
    deriv1 <- function(x, parm){
      d1 <- rep(1, length(x) )
      d2 <- x
      cbind(d1, d2)
   }
    ## Defining the ED function

    ## Defining the inverse function

    ## Defining descriptive text
    text <- "Straight line"

    ## Returning the function with self starter and names
    returnList <- list(fct = fct, ssfct = ssfct, deriv1 = deriv1, names = pnames, text = text, noParm = sum(is.na(fixed)))

    class(returnList) <- "drcMean"
    invisible(returnList)
}

#linear model with no intercept#########################################
linearOrigin.fun <- function(predictor, b) {
                      b * predictor
}

linearOrigin.Init <- function(mCall, LHS, data) {
          xy <- sortedXyData(mCall[["predictor"]], LHS, data)
          lmFit <- lm((xy[, "y"]) ~ xy[, "x"]-1)
          coefs <- coef(lmFit)
          b <- coefs[1]
          value <- c(b)
          names(value) <- mCall[c("b")]
          value
}

NLS.linearOrigin <- selfStart(linearOrigin.fun, linearOrigin.Init, parameters=c("b"))

"DRC.linearOrigin" <- function(fixed = c(NA), names = c("b"))
{
    ## Checking arguments
    numParm <- 1
    if (!is.character(names) | !(length(names) == numParm))
      {stop("Not correct 'names' argument")}
    if (!(length(fixed) == numParm))
      {stop("Not correct 'fixed' argument")}

    ## Fixing parameters (using argument 'fixed')
    notFixed <- is.na(fixed)
    parmVec <- rep(0, numParm)
    parmVec[!notFixed] <- fixed[!notFixed]

    ## Defining the non-linear function
    fct <- function(x, parm)
    {
        parmMat <- matrix(parmVec, nrow(parm), numParm, byrow = TRUE)
        parmMat[, notFixed] <- parm

        b <- parmMat[, 1]
        linearOrigin.fun(x, b)
    }

    ## Defining self starter function
    ssfct <- function(dataf)
    {
        x <- dataf[, 1]
        y <- dataf[, 2]

        #regression on pseudo y values
        pseudoY <- y
        pseudoX <- x
        coefs <- coef( lm(pseudoY ~ pseudoX - 1) )
        b <- coefs[1]

        return(c(b)[notFixed])
    }

    ## Defining names
    pnames <- names[notFixed]

    ## Defining derivatives
    deriv1 <- function(x, parm){
      d1 <- x
      cbind(d1)
   }
    ## Defining the ED function

    ## Defining the inverse function

    ## Defining descriptive text
    text <- "Straight line through origin"

    ## Returning the function with self starter and names
    returnList <- list(fct = fct, ssfct = ssfct, deriv1 = deriv1, names = pnames, text = text, noParm = sum(is.na(fixed)))

    class(returnList) <- "drcMean"
    invisible(returnList)
}

#Polynomial regression (2nd order) ##############################################
poly2.fun <- function(predictor, b) {
                      a + b * predictor + c * (predictor^2)
}

poly2.Init <- function(mCall, LHS, data) {
          xy <- sortedXyData(mCall[["predictor"]], LHS, data)
          lmFit <- lm((xy[, "y"]) ~ xy[, "x"] + I( xy[, "x"]^2))
          coefs <- coef(lmFit)
          a <- coefs[1]
          b <- coefs[2]
          c <- coefs[3]
          value <- c(a, b, c)
          names(value) <- mCall[c("a", "b", "c")]
          value
}
NLS.poly2 <- selfStart(poly2.fun, poly2.Init, parameters=c("a", "b", "c"))

"DRC.poly2" <- function(fixed = c(NA, NA, NA), names = c("a", "b", "c"))
{
  ## Checking arguments
  numParm <- 3
  if (!is.character(names) | !(length(names) == numParm)) {stop("Not correct 'names' argument")}
  if (!(length(fixed) == numParm)) {stop("Not correct 'fixed' argument")}

  ## Fixing parameters (using argument 'fixed')
  notFixed <- is.na(fixed)
  parmVec <- rep(0, numParm)
  parmVec[!notFixed] <- fixed[!notFixed]

  ## Defining the non-linear function
  fct <- function(x, parm)
  {
    parmMat <- matrix(parmVec, nrow(parm), numParm, byrow = TRUE)
    parmMat[, notFixed] <- parm

    a <- parmMat[, 1]; b <- parmMat[, 2]; c <- parmMat[, 3]
    a + b * x + c * (x^2)
  }

  ## Defining self starter function
  ssfct <- function(dataf)
  {
    x <- dataf[, 1]
    y <- dataf[, 2]

    #regression on pseudo y values
    pseudoY <- y
    pseudoX <- x
    coefs <- coef( lm(pseudoY ~ pseudoX + I(pseudoX^2)) )
    a <- coefs[1]
    b <- coefs[2]
    c <- coefs[3]

    return(c(a, b, c)[notFixed])
  }

  ## Defining names
  pnames <- names[notFixed]

  ## Defining derivatives
  deriv1 <- function(x, parm){
      d1 <- rep(1, length(x) )
      d2 <- x
      d3 <- x^2
      cbind(d1, d2, d3)
   }
  ## Defining the ED function

  ## Defining the inverse function

  ## Defining descriptive text
  text <- "Second Order Polynomial"

  ## Returning the function with self starter and names
  returnList <- list(fct = fct, ssfct = ssfct, names = pnames, 
    text = text, noParm = sum(is.na(fixed)), deriv1 = deriv1)

  class(returnList) <- "drcMean"
  invisible(returnList)
}

# Exponential growth ####################################################
EXD.fun <- function(X, c, d, e){
  c + (d - c) * exp( X / e)  
}

expoGrowth.fun <- function(predictor, a, b) {
                      a * exp(b * predictor)
}

expoGrowth.Init <- function(mCall, LHS, data) {
          xy <- sortedXyData(mCall[["predictor"]], LHS, data)
          x <-  xy[, "x"]; y <- xy[, "y"]
          ## Linear regression on pseudo y values
          pseudoY <- log(y + 0.000001)
          coefs <- coef( lm(pseudoY ~ x) )
          k <- coefs[1]
          b <- coefs[2]
          a <- exp(k)
          value <- c(a, b)
          names(value) <- mCall[c("init", "k")]
          value
}

NLS.expoGrowth <- selfStart(expoGrowth.fun, expoGrowth.Init, parameters=c("init", "k"))

"DRC.expoGrowth" <-
function(fixed = c(NA, NA), names = c("init", "k"))
{
    ## Checking arguments
    numParm <- 2
    if (!is.character(names) | !(length(names) == numParm)) {stop("Not correct 'names' argument")}
    if (!(length(fixed) == numParm)) {stop("Not correct 'fixed' argument")}

    ## Fixing parameters (using argument 'fixed')
    notFixed <- is.na(fixed)
    parmVec <- rep(0, numParm)
    parmVec[!notFixed] <- fixed[!notFixed]

    ## Defining the non-linear function
    fct <- function(x, parm)
    {
        parmMat <- matrix(parmVec, nrow(parm), numParm, byrow = TRUE)
        parmMat[, notFixed] <- parm

        init <- parmMat[, 1]; k <- parmMat[, 2]
        init * exp ( k * x )
    }

    ## Defining self starter function
    ssfct <- function(dataf)
    {
        x <- dataf[, 1]
        y <- dataf[, 2]

        ## Linear regression on pseudo y values
        pseudoY <- log( y + 0.000001)
        coefs <- coef( lm(pseudoY ~ x) )
        init <- exp(coefs[1])
        m <- coefs[2]

        return(c(init, m)[notFixed])
    }

    ## Defining names
    pnames <- names[notFixed]

    ## Defining derivatives

    ## Defining the ED function

    ## Defining the inverse function

    ## Defining descriptive text
    text <- "Exponential Growth Model"

    ## Returning the function with self starter and names
    returnList <- list(fct = fct, ssfct = ssfct, names = pnames, text = text, noParm = sum(is.na(fixed)))

    class(returnList) <- "drcMean"
    invisible(returnList)
}



#Negative exponential ###########################################################
negExpDist.fun <- function(predictor, c) {
                      x <- predictor
                      1 - exp (- c * x)
}

negExpDist.Init <- function(mCall, LHS, data) {
          xy <- sortedXyData(mCall[["predictor"]], LHS, data)
          x <-  xy[, "x"]; y <- xy[, "y"]
          ## Linear regression on pseudo y values
          pseudoY <- log( 1 - y )
          coefs <- coef( lm(pseudoY ~ x - 1) )
          c <- - coefs[1]
          value <- c(c)
          names(value) <- mCall[c("c")]
          value
}

NLS.negExpDist <- selfStart(negExpDist.fun, negExpDist.Init, parameters=c("c"))

"DRC.negExpDist" <-
function(fixed = c(NA, NA), names = c("c"))
{
    ## Checking arguments
    numParm <- 1
    if (!is.character(names) | !(length(names) == numParm)) {stop("Not correct 'names' argument")}
    if (!(length(fixed) == numParm)) {stop("Not correct 'fixed' argument")}

    ## Fixing parameters (using argument 'fixed')
    notFixed <- is.na(fixed)
    parmVec <- rep(0, numParm)
    parmVec[!notFixed] <- fixed[!notFixed]

    ## Defining the non-linear function
    fct <- function(x, parm)
    {
        parmMat <- matrix(parmVec, nrow(parm), numParm, byrow = TRUE)
        parmMat[, notFixed] <- parm

        c <- parmMat[, 1]
        1 - exp (- c * x)
    }

    ## Defining self starter function
    ssfct <- function(dataf)
    {
        x <- dataf[, 1]
        y <- dataf[, 2]

         ## Linear regression on pseudo y values
        pseudoY <- log( 1 - y )
        coefs <- coef( lm(pseudoY ~ x - 1) )
        c <- - coefs[1]

        return(c(c)[notFixed])
    }

    ## Defining names
    pnames <- names[notFixed]

    ## Defining derivatives

    ## Defining the ED function

    ## Defining the inverse function

    ## Defining descriptive text
    text <- "Negative exponential cumulative distribution function"

    ## Returning the function with self starter and names
    returnList <- list(fct = fct, ssfct = ssfct, names = pnames, text = text, noParm = sum(is.na(fixed)))

    class(returnList) <- "drcMean"
    invisible(returnList)
}


#Power Curve ########################################################
powerCurve.fun <- function(predictor, a, b) {
                      a * ( predictor ^ b )
}

powerCurve.Init <- function(mCall, LHS, data) {
          xy <- sortedXyData(mCall[["predictor"]], LHS, data)
          lmFit <- lm(log(xy[, "y"]) ~ log(xy[, "x"]))
          coefs <- coef(lmFit)
          a <- exp(coefs[1])
          b <- coefs[2]
          value <- c(a, b)
          names(value) <- mCall[c("a", "b")]
          value
}

NLS.powerCurve <- selfStart(powerCurve.fun, powerCurve.Init, parameters=c("a", "b"))

"DRC.powerCurve" <- function(fixed = c(NA, NA), names = c("a", "b"))
{
    ## Checking arguments
    numParm <- 2
    if (!is.character(names) | !(length(names) == numParm)) {stop("Not correct 'names' argument")}
    if (!(length(fixed) == numParm)) {stop("Not correct 'fixed' argument")}

    ## Fixing parameters (using argument 'fixed')
    notFixed <- is.na(fixed)
    parmVec <- rep(0, numParm)
    parmVec[!notFixed] <- fixed[!notFixed]

    ## Defining the non-linear function
    fct <- function(x, parm)
    {


        parmMat <- matrix(parmVec, nrow(parm), numParm, byrow = TRUE)
        parmMat[, notFixed] <- parm

        a <- parmMat[, 1]; b <- parmMat[, 2]
        a * x ^(b)
    }

    ## Defining self starter function
    ssfct <- function(dataf)
    {
        x <- dataf[, 1]
        y <- dataf[, 2]

        #regression on pseudo y values
        pseudoY <- log( y + 0.00001)
        pseudoX <- log(x)
        coefs <- coef( lm(pseudoY ~ pseudoX) )
        a <- exp(coefs[1])

        b <- coefs[2]

        return(c(a, b)[notFixed])
    }

    ## Defining names
    pnames <- names[notFixed]

    ## Defining derivatives

    ## Defining the ED function

    ## Defining the inverse function

    ## Defining descriptive text
    text <- "Power curve (Freundlich equation)"

    ## Returning the function with self starter and names
    returnList <- list(fct = fct, ssfct = ssfct, names = pnames, text = text, noParm = sum(is.na(fixed)))

    class(returnList) <- "drcMean"
    invisible(returnList)
}

#Logarithmic regression ########################################
logCurve.fun <- function(predictor, a, b) {
                      x <- predictor
                      a  + b * log(x)
}

logCurve.Init <- function(mCall, LHS, data) {
          xy <- sortedXyData(mCall[["predictor"]], LHS, data)
          x <-  xy[, "x"]; y <- xy[, "y"]
          pseudoY <-  y 
          pseudoX <- log(x)
          coefs <- coef( lm(pseudoY ~ pseudoX) )
          a <- coefs[1]        
          b <- coefs[2]
          value <- c(a, b)
          names(value) <- mCall[c("a", "b")]
          value
}

NLS.logCurve <- selfStart(logCurve.fun, logCurve.Init, parameters=c("a", "b"))

"DRC.logCurve" <-
function(fixed = c(NA, NA), names = c("a", "b"))
{
    ## Checking arguments
    numParm <- 2
    if (!is.character(names) | !(length(names) == numParm)) {stop("Not correct 'names' argument")}
    if (!(length(fixed) == numParm)) {stop("Not correct 'fixed' argument")}

    ## Fixing parameters (using argument 'fixed')
    notFixed <- is.na(fixed)
    parmVec <- rep(0, numParm)
    parmVec[!notFixed] <- fixed[!notFixed]

    ## Defining the non-linear function
    fct <- function(x, parm)
    {


        parmMat <- matrix(parmVec, nrow(parm), numParm, byrow = TRUE)
        parmMat[, notFixed] <- parm

        a <- parmMat[, 1]; b <- parmMat[, 2]
        a  + b * log(x)
    }

    ## Defining self starter function
    ssfct <- function(dataf)
    {
        x <- dataf[, 1]
        y <- dataf[, 2]

        #regression on pseudo y values
        pseudoY <-  y
        pseudoX <- log(x)
        coefs <- coef( lm(pseudoY ~ pseudoX) )
        a <- coefs[1]
        b <- coefs[2]

        return(c(a, b)[notFixed])
    }

    ## Defining names
    pnames <- names[notFixed]

    ## Defining derivatives

    ## Defining the ED function

    ## Defining the inverse function

    ## Defining descriptive text
    text <- "Linear regression on log-transformed x"

    ## Returning the function with self starter and names
    returnList <- list(fct = fct, ssfct = ssfct, names = pnames, text = text, noParm = sum(is.na(fixed)))

    class(returnList) <- "drcMean"
    invisible(returnList)
}

#Logarithmic regression without intercept ################################
logCurveNI.fun <- function(predictor, b) {
                      x <- predictor
                      b * log(x)
}

logCurveNI.Init <- function(mCall, LHS, data) {
          xy <- sortedXyData(mCall[["predictor"]], LHS, data)
          x <-  xy[, "x"]; y <- xy[, "y"]
          pseudoY <-  y 
          pseudoX <- log(x)
          coefs <- coef( lm(pseudoY ~ pseudoX - 1) )
          b <- coefs[1]
          value <- c(b)
          names(value) <- mCall[c("b")]
          value
}

NLS.logCurveNI <- selfStart(logCurveNI.fun, logCurveNI.Init, parameters=c("b"))

#Yield-Loss function (rectangular hyperbola) ######################
YL.fun <- function(predictor, i, A) {
            i * predictor/(1 + i/A * predictor)
}

YL.Init <- function(mCall, LHS, data) {
          xy <- sortedXyData(mCall[["predictor"]], LHS, data)
          x <- xy[, "x"]; y <- xy[, "y"]
          pseudoX <- 1 / x[x > 0]; pseudoY <- 1 / y[x > 0]
          lmFit <- lm(pseudoY ~ pseudoX)
          coefs <- coef(lmFit)
          A <- 1 / coefs[1]
          i <- 1 / coefs[2]
          value <- c(i, A)
          names(value) <- mCall[c("i", "A")]
          value
}
NLS.YL <- selfStart(YL.fun, YL.Init, parameters=c("i", "A"))

"DRC.YL" <- function(fixed = c(NA, NA), names = c("i", "A")) {
    ## Checking arguments
    numParm <- 2
    if (!is.character(names) | !(length(names) == numParm)) {stop("Not correct 'names' argument")}
    if (!(length(fixed) == numParm)) {stop("Not correct 'fixed' argument")}

    ## Fixing parameters (using argument 'fixed')
    notFixed <- is.na(fixed)
    parmVec <- rep(0, numParm)
    parmVec[!notFixed] <- fixed[!notFixed]

    ## Defining the non-linear function
    fct <- function(x, parm) {
      parmMat <- matrix(parmVec, nrow(parm), numParm, byrow = TRUE)
      parmMat[, notFixed] <- parm
      i <- parmMat[, 1]; A <- parmMat[, 2]
      YL.fun(x, i, A)
    }

    ## Defining self starter function
    ssfct <- function(dataf) {
      x <- dataf[, 1]
      y <- dataf[, 2]

      #regression on pseudo y values
      pseudoY <- 1 /  y[x > 0] 
      pseudoX <- 1 / x [x > 0] 
      coefs <- coef( lm(pseudoY ~ pseudoX) )
      A <- 1 / coefs[1]; i <- 1 / coefs[2]

      return(c(i, A)[notFixed])
    }

    ## Defining names
    pnames <- names[notFixed]

    ## Defining derivatives

    ## Defining the ED function

    ## Defining the inverse function

    ## Defining descriptive text
    text <- "Yield-Loss function (Cousens, 1985)"

    ## Returning the function with self starter and names
    returnList <- list(fct = fct, ssfct = ssfct, names = pnames, text = text, noParm = sum(is.na(fixed)))

    class(returnList) <- "drcMean"
    invisible(returnList)
}

#Yield - Weed Density function ###################
cousens85.fun <- function(Dens, Ywf, i, A) {
      Ywf*(1 - (i * Dens) / (100 * (1 + i * Dens/A)))
}

"DRC.cousens85" <-
  function(fixed = c(NA, NA, NA), names = c("YWF", "i", "a"))
  {
    ## Checking arguments
    numParm <- 3
    if (!is.character(names) | !(length(names) == numParm)) {stop("Not correct 'names' argument")}
    if (!(length(fixed) == numParm)) {stop("Not correct 'fixed' argument")}

    ## Fixing parameters (using argument 'fixed')
    notFixed <- is.na(fixed)
    parmVec <- rep(0, numParm)
    parmVec[!notFixed] <- fixed[!notFixed]

    ## Defining the non-linear function
    fct <- function(x, parm)
    {
      parmMat <- matrix(parmVec, nrow(parm), numParm, byrow = TRUE)
      parmMat[, notFixed] <- parm

      YWF <- parmMat[, 1]; i <- parmMat[, 2]; a <- parmMat[, 3]
      YWF*(1 - (i*x) / (100 * (1 + i * x/a)))
    }

    ## Defining self starter function
    ssfct <- function(dataf)
    {
      x <- dataf[, 1]
      y <- dataf[, 2]

      YWF <- max(y)+10e-06
      YL <- (1 - y/YWF)*100
      #regression on pseudo y values
      pseudoY <- 1 /  (YL + 0.000001)
      pseudoX <- 1 / (x + 0.00001)
      coefs <- coef( lm(pseudoY ~ pseudoX) )
      a <- 1 / coefs[1]; i <- 1 / coefs[2]

      return(c(YWF, i, a)[notFixed])
    }

    ## Defining names
    pnames <- names[notFixed]

    ## Defining derivatives

    ## Defining the ED function

    ## Defining the inverse function

    ## Defining descriptive text
    text <- "Yield-Weed Density function (Cousens, 1985)"

    ## Returning the function with self starter and names
    returnList <- list(fct = fct, ssfct = ssfct, names = pnames, text = text, noParm = sum(is.na(fixed)))

    class(returnList) <- "drcMean"
    invisible(returnList)
}

#Logistic curve with two-parameters
"L.2" <-
function(fixed = c(NA, NA), names = c("b", "e"), ...)
{
    ## Checking arguments
    numParm <- 2
    
    if (!is.character(names) | !(length(names)==numParm)) {stop("Not correct names argument")}
    if (!(length(fixed)==numParm)) {stop("Not correct length of 'fixed' argument")}

    return(logistic(fixed = c(fixed[1], 0, 1, fixed[2], 1), names = c(names[1], "c", "d", names[2], "f"), 
    fctName = as.character(match.call()[[1]]), 
    fctText = "Logistic (ED50 as parameter) with lower limit fixed at 0 and upper limit fixed to 1", ...))
}

L.fun <- function(X, b, c, d, e) {
  c + (d - c)/(1 + exp ( b * ( X - e)))
}

G.fun <- function(X, b, c, d, e) {
  c + (d - c) * (exp (- exp ( b * ( X - e))))
}

E.fun <- function(X, b, c, d, e) {
  c + (d - c) * (1 - exp( - exp ( - b * ( X - e))))
}

#Rational function ################################################
Rational.fun <- function(predictor, a, b, c) {
                      x <- predictor
                      (b + c*x) / (1 + a*x)
}

Rational.Init <- function(mCall, LHS, data) {
          xy <- sortedXyData(mCall[["predictor"]], LHS, data)
          x <-  xy[, "x"]; y <- xy[, "y"]
          pseudoY <-  y 
          pseudoX <- x
          pseudoXY <- x*y
          coefs <- coef( lm(pseudoY ~ pseudoX + pseudoXY) )
          b <- coefs[1]        
          c <- coefs[2]
          a <- - coefs[3]
          value <- c(a, b, c)
          names(value) <- mCall[c("a", "b", "c")]
          value
}

NLS.Rational <- selfStart(Rational.fun, Rational.Init, parameters=c("a", "b", "c"))

"DRC.Rational" <-
function(fixed = c(NA, NA, NA), names = c("a", "b", "c"))
{
    ## Checking arguments
    numParm <- 3
    if (!is.character(names) | !(length(names) == numParm)) {stop("Not correct 'names' argument")}
    if (!(length(fixed) == numParm)) {stop("Not correct 'fixed' argument")}

    ## Fixing parameters (using argument 'fixed')
    notFixed <- is.na(fixed)
    parmVec <- rep(0, numParm)
    parmVec[!notFixed] <- fixed[!notFixed]

    ## Defining the non-linear function
    fct <- function(x, parm)
    {
        parmMat <- matrix(parmVec, nrow(parm), numParm, byrow = TRUE)
        parmMat[, notFixed] <- parm

        a <- parmMat[, 1]; b <- parmMat[, 2]; c <- parmMat[, 3]
       (b + c*x)/(1 + a*x)
    }

    ## Defining self starter function
    ssfct <- function(dataf)
    {
        x <- dataf[, 1]
        y <- dataf[, 2]
        xy <- dataf[,1]*dataf[,2]

        ## Linear regression on pseudo y values

        coefs <- coef(lm(y ~ x + xy) )
        b <- coefs[1]
        c <- coefs[2]
        a <- - coefs[3]

        return(c(a, b, c)[notFixed])
    }

    ## Defining names
    pnames <- names[notFixed]

    ## Defining derivatives

    ## Defining the ED function

    ## Defining the inverse function

    ## Defining descriptive text
    text <- "Inverse polynomial 1"

    ## Returning the function with self starter and names
    returnList <- list(fct = fct, ssfct = ssfct, names = pnames, text = text, noParm = sum(is.na(fixed)))

    class(returnList) <- "drcMean"
    invisible(returnList)
}





