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




# Exponential growth ####################################################
EXD.fun <- function(X, c, d, e){
  c + (d - c) * exp( X / e)  
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


#Rational function ################################################
# Ratio of two polynomials ############################
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





