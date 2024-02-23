"DRClogiGrowth.1" <-
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
    
        a <- parmMat[, 1]; b <- parmMat[, 3]; c <- parmMat[, 2] 
        a / (1 + exp(- b * x + c))
    }

    ## Defining self starter function        
    ssfct <- function(dataf) 
    {       
        x <- dataf[, 1]
        y <- dataf[, 2]

        a <- max(y) * 1.05
        
        ## Linear regression on pseudo y values
        pseudoY <- log( (a / y ) - 1 )
        coefs <- coef( lm(pseudoY ~ x) )
        b <- - coefs[2]
        c <- coefs[1]
	    
        return(c(a, b, c)[notFixed])
    }

    ## Defining names
    pnames <- names[notFixed]

    ## Defining derivatives

    ## Defining the ED function
       
    ## Defining the inverse function
    
    ## Defining descriptive text
    text <- "Logistic Growth Model - 1"    

    ## Returning the function with self starter and names
    returnList <- list(fct = fct, ssfct = ssfct, names = pnames, text = text, noParm = sum(is.na(fixed))) 
    
    class(returnList) <- "drcMean"
    invisible(returnList)
}

"DRClogiGrowth.2" <-
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
    
        a <- parmMat[, 1]; b <- parmMat[, 3]; c <- parmMat[, 2] 
        a / (1 + b * exp(- c * x))
    }

    ## Defining self starter function        
    ssfct <- function(dataf) 
    {       
        x <- dataf[, 1]
        y <- dataf[, 2]

        a <- max(y) * 1.05
        
        ## Linear regression on pseudo y values
        pseudoY <- log( (a / y ) - 1 )
        coefs <- coef( lm(pseudoY ~ x) )
        c <- - coefs[2]
        b <- exp(coefs[1])
	    
        return(c(a, b, c)[notFixed])
    }

    ## Defining names
    pnames <- names[notFixed]

    ## Defining derivatives

    ## Defining the ED function
       
    ## Defining the inverse function
    
    ## Defining descriptive text
    text <- "Logistic Growth Model - 2"    

    ## Returning the function with self starter and names
    returnList <- list(fct = fct, ssfct = ssfct, names = pnames, text = text, noParm = sum(is.na(fixed))) 
    
    class(returnList) <- "drcMean"
    invisible(returnList)
}

"DRClogiGrowth.3" <-
function(fixed = c(NA, NA, NA), names = c("init", "m", "plateau"))
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
    
        W0 <- parmMat[, 1]; Wf <- parmMat[, 3]; m <- parmMat[, 2] 
        W0 * Wf / (W0 + (Wf - W0) * exp( - m * x)) 
    }

    ## Defining self starter function        
    ssfct <- function(dataf) 
    {       
        x <- dataf[, 1]
        y <- dataf[, 2]

        plateau <- max(y) * 1.05
        
        ## Linear regression on pseudo y values
        pseudoY <- log( (plateau / (y + 0.0001) ) - 1 )
        coefs <- coef( lm(pseudoY ~ x) )
        b <- exp(coefs[1])
        init <- plateau / (1 + b)
        m <- - coefs[2]
	    
        return(c(init, m, plateau)[notFixed])
    }

    ## Defining names
    pnames <- names[notFixed]

    ## Defining derivatives

    ## Defining the ED function
       
    ## Defining the inverse function
    
    ## Defining descriptive text
    text <- "Logistic Growth Model - 3"    

    ## Returning the function with self starter and names
    returnList <- list(fct = fct, ssfct = ssfct, names = pnames, text = text, noParm = sum(is.na(fixed))) 
    
    class(returnList) <- "drcMean"
    invisible(returnList)
}

"DRClogiGrowth.4" <-
function(fixed = c(NA, NA, NA), names = c("m", "plateau", "t50"))
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
            
        m <- parmMat[, 1]; plateau <- parmMat[, 2]; t50 <- parmMat[, 3] 
        plateau / (1 + exp(- m * (x - t50))) 
    }

    ## Defining self starter function        
    ssfct <- function(dataf) 
    {       
        x <- dataf[, 1]
        y <- dataf[, 2]

        plateau <- max(y) * 1.05
        
        ## Linear regression on pseudo y values
        pseudoY <- log( (plateau / (y + 0.0001) ) - 1 )
        coefs <- coef( lm(pseudoY ~ x) )
        m <- - coefs[2]
        t50 <- coefs[1] / m
	    
        return(c(m, plateau, t50)[notFixed])
    }

    ## Defining names
    pnames <- names[notFixed]

    ## Defining derivatives

    ## Defining the ED function
       
    ## Defining the inverse function
    
    ## Defining descriptive text
    text <- "Logistic Growth Model - 4"    

    ## Returning the function with self starter and names
    returnList <- list(fct = fct, ssfct = ssfct, names = pnames, text = text, noParm = sum(is.na(fixed))) 
    
    class(returnList) <- "drcMean"
    invisible(returnList)
}

"DRCfpl" <-
function(fixed = c(NA, NA, NA, NA), names = c("c", "d", "mu", "sigma"))
{
    ## Checking arguments
    numParm <- 4
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
            
        c <- parmMat[, 1]; d <- parmMat[, 2]; mu <- parmMat[, 3]; sigma <- parmMat[, 3] 
        c + (d - c) / (1 + exp((mu-x)/sigma)) 
    }

    ## Defining self starter function        
    ssfct <- function(dataf) 
    {       
        x <- dataf[, 1]
        y <- dataf[, 2]

        d <- max(y) * 1.05
        c <- min(y) - 0.0001
        ## Linear regression on pseudo y values
        pseudoY <- log((d - y) / (y - c ))
        coefs <- coef( lm(pseudoY ~ x) )
        m <- - coefs[2]
        mu <- coefs[1] / m
        sigma <- 1/m
	    
        return(c(c, d, mu, sigma)[notFixed])
    }

    ## Defining names
    pnames <- names[notFixed]

    ## Defining derivatives

    ## Defining the ED function
       
    ## Defining the inverse function
    
    ## Defining descriptive text
    text <- "Logistic function (4 parms)"    

    ## Returning the function with self starter and names
    returnList <- list(fct = fct, ssfct = ssfct, names = pnames, text = text, noParm = sum(is.na(fixed))) 
    
    class(returnList) <- "drcMean"
    invisible(returnList)
}

"logLogisticD.2" <-
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
        1 - (1 / (1 + exp(a + b*log(x)))) 
    }

    ## Defining self starter function        
    ssfct <- function(dataf) 
    {       
        x <- dataf[, 1]
        y <- 1 - dataf[, 2] - 0.000001

        ## Linear regression on pseudo y values
        pseudoY <- log( y / (1-y ))
        coefs <- coef( lm(pseudoY ~ x) )
        b <- coefs[2]
        a <- coefs[1]
	    
        return(c(a, b)[notFixed])
    }

    ## Defining names
    pnames <- names[notFixed]

    ## Defining derivatives

    ## Defining the ED function
       
    ## Defining the inverse function
    
    ## Defining descriptive text
    text <- "Traditional log-Logistic model"    

    ## Returning the function with self starter and names
    returnList <- list(fct = fct, ssfct = ssfct, names = pnames, text = text, noParm = sum(is.na(fixed))) 
    
    class(returnList) <- "drcMean"
    invisible(returnList)
}

"logLogisticD.3" <-
function(fixed = c(NA, NA, NA), names = c("plateau","a", "b"))
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
            
        plateau <- parmMat[,1]; a <- parmMat[, 2]; b <- parmMat[, 3]
        1 - (plateau / (1 + exp(a + b*log(x)))) 
    }

    ## Defining self starter function        
    ssfct <- function(dataf) 
    {       
        x <- dataf[, 1]
        y <- 1 - dataf[, 2] + 0.000001

        ## Linear regression on pseudo y values
        plateau <- max(y) * 1.05
        
        pseudoY <- log( y / (1-y))
        coefs <- coef( lm(pseudoY ~ x) )
        b <- coefs[2]
        a <- coefs[1]
	    
        return(c(plateau, a, b)[notFixed])
    }

    ## Defining names
    pnames <- names[notFixed]

    ## Defining derivatives

    ## Defining the ED function
       
    ## Defining the inverse function
    
    ## Defining descriptive text
    text <- "Traditional log-Logistic model with upper asymptote"    

    ## Returning the function with self starter and names
    returnList <- list(fct = fct, ssfct = ssfct, names = pnames, text = text, noParm = sum(is.na(fixed))) 
    
    class(returnList) <- "drcMean"
    invisible(returnList)
}

"logLogistic.2" <-
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
        1 / (1 + exp(-(a + b*log(x)))) 
    }

    ## Defining self starter function        
    ssfct <- function(dataf) 
    {       
        x <- dataf[, 1]
        y <- dataf[, 2] + 0.000001

        ## Linear regression on pseudo y values
        pseudoY <- log( y / (1-y) )
        coefs <- coef( lm(pseudoY ~ x) )
        b <- coefs[2]
        a <- coefs[1]
	    
        return(c(a, b)[notFixed])
    }

    ## Defining names
    pnames <- names[notFixed]

    ## Defining derivatives

    ## Defining the ED function
       
    ## Defining the inverse function
    
    ## Defining descriptive text
    text <- "Traditional log-Logistic model"    

    ## Returning the function with self starter and names
    returnList <- list(fct = fct, ssfct = ssfct, names = pnames, text = text, noParm = sum(is.na(fixed))) 
    
    class(returnList) <- "drcMean"
    invisible(returnList)
}

"logLogistic.3" <-
function(fixed = c(NA, NA, NA), names = c("plateau","a", "b"))
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
            
        plateau <- parmMat[,1]; a <- parmMat[, 2]; b <- parmMat[, 3]
        plateau / (1 + exp(-(a + b*log(x)))) 
    }

    ## Defining self starter function        
    ssfct <- function(dataf) 
    {       
        x <- dataf[, 1]
        y <- dataf[, 2] + 0.000001

        ## Linear regression on pseudo y values
        plateau <- max(y) * 1.05
        
        pseudoY <- log( y / (plateau - y) )
        coefs <- coef( lm(pseudoY ~ x) )
        b <- coefs[2]
        a <- coefs[1]
	    
        return(c(plateau, a, b)[notFixed])
    }

    ## Defining names
    pnames <- names[notFixed]

    ## Defining derivatives

    ## Defining the ED function
       
    ## Defining the inverse function
    
    ## Defining descriptive text
    text <- "Traditional log-Logistic model with upper asymptote"    

    ## Returning the function with self starter and names
    returnList <- list(fct = fct, ssfct = ssfct, names = pnames, text = text, noParm = sum(is.na(fixed))) 
    
    class(returnList) <- "drcMean"
    invisible(returnList)
}
