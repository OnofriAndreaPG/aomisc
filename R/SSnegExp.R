#Negative exponential Function #################################################
negExp.fun <- function(predictor, a, c) {
                      x <- predictor
                      a * (1 - exp (- c * x))
}

negExp.Init <- function(mCall, LHS, data) {
          xy <- sortedXyData(mCall[["predictor"]], LHS, data)
          x <-  xy[, "x"]; y <- xy[, "y"]
          plateau <- max(y) * 1.05        
          ## Linear regression on pseudo y values
          pseudoY <- log( 1 - (y / plateau ) )
          coefs <- coef( lm(pseudoY ~ x - 1) )
          a <- plateau
          c <- - coefs[1]
          value <- c(a, c)
          names(value) <- mCall[c("a", "c")]
          value
}

NLS.negExp <- selfStart(negExp.fun, negExp.Init, parameters=c("a", "c"))

"DRC.negExp" <-
function(fixed = c(NA, NA), names = c("a", "c"))
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

        a <- parmMat[, 1]; c <- parmMat[, 2]
        a * (1 - exp (- c * x))
    }

    ## Defining self starter function
    ssfct <- function(dataf)
    {
        x <- dataf[, 1]
        y <- dataf[, 2]

        plateau <- max(y) * 1.05

        ## Linear regression on pseudo y values
        pseudoY <- log( 1 - (y / plateau ) )
        coefs <- coef( lm(pseudoY ~ x - 1) )
        a <- plateau
        c <- - coefs[1]

        return(c(a, c)[notFixed])
    }

    ## Defining names
    pnames <- names[notFixed]

    ## Defining derivatives

    ## Defining the ED function

    ## Defining the inverse function

    ## Defining descriptive text
    text <- "Negative exponential function"

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



