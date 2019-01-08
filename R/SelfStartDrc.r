"linear" <- function(fixed = c(NA, NA), names = c("a", "b"))
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
        a + b * x
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
    returnList <- list(fct = fct, deriv1= deriv1, ssfct = ssfct, names = pnames, text = text, noParm = sum(is.na(fixed)))

    class(returnList) <- "drcMean"
    invisible(returnList)
}

"parabolic" <- function(fixed = c(NA, NA, NA), names = c("a", "b", "c"))
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

  ## Defining the ED function

  ## Defining the inverse function

  ## Defining descriptive text
  text <- "Second Order Polynomial"

  ## Returning the function with self starter and names
  returnList <- list(fct = fct, ssfct = ssfct, names = pnames, text = text, noParm = sum(is.na(fixed)))

  class(returnList) <- "drcMean"
  invisible(returnList)
}

"firstOrder" <-
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
        init * exp ( - k * x )
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
        m <- - coefs[2]

        return(c(init, m)[notFixed])
    }

    ## Defining names
    pnames <- names[notFixed]

    ## Defining derivatives

    ## Defining the ED function

    ## Defining the inverse function

    ## Defining descriptive text
    text <- "Exponential Decay Model"

    ## Returning the function with self starter and names
    returnList <- list(fct = fct, ssfct = ssfct, names = pnames, text = text, noParm = sum(is.na(fixed)))

    class(returnList) <- "drcMean"
    invisible(returnList)
}



