
# Exponential growth ####################################################
EXD.fun <- function(X, c, d, e){
  c + (d - c) * exp( X / e)  
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





