"boxcox.nls" <- function(object, lambda = seq(-2, 2, 1/10), plotit = TRUE, start,
eps = 1/50, bcAdd = 0, level = 0.95, xlab = expression(lambda), ylab = "log-likelihood", ...)
{
    ## Defining the Box-Cox modified power transformations
    bfFct <- function(lv)
    {
        function(x) {if (abs(lv) < eps) {return(log(x+bcAdd))} else {return(((x+bcAdd)^lv-1)/lv)}}    
    }

    evalForm <- eval(formula(object)[[3]][[1]])
    bfFct2 <- function(lv)
    {
        bcFct2 <- function(x)
        {
            ## Transforming the mean
            bcFct <- bfFct(lv)
            bcFctVal <- bcFct(x)

#            ## Assigning sttributes
#            attrList <- attributes(evalForm)
#            lenAL <- length(attrList)
#            namAL <- names(attrList)
#     
#            for (i in 1:lenAL)
#            {
##                attr(bcFctVal, namAL[i]) <- attrList[[i]]
#            }
#            
            ## Adjusting gradient
            grad1 <- attr(bcFctVal, "gradient")
            grad2 <- grad1*(x^(lv-1))
            attr(bcFctVal, "gradient") <- grad2 
            
            bcFctVal
        }
        bcFct2
     }
     
     if (!inherits(evalForm , "selfStart"))
     {
         bfFct2 <- bfFct
     }    
         
#    bcFct <- function(x) {if (abs(lv) < eps) {return(log(x+bcAdd))} else {return(((x+bcAdd)^lv-1)/lv)}}
#    assign("bcFct", bcFct, envir = .GlobalEnv)
#    newFormula <- bcFct(.) ~ bcFct(.)
    if (missing(start)) 
    {
        startVec <- coef(object)
    } else {
        startVec <- start
    }
    
    ## Defining the log-likelihood function
    llFct <- function(object, lv)
    {
        sumObj <- summary(object)
        N <- sum(sumObj$df)
        Ji <- (eval(object$data)[, as.character(formula(object)[[2]])[2]])^(lv-1)
        
        -N*log(sqrt(sum(residuals(object)^2)/N))-N/2+sum(log(Ji))
    }
    
    lenlam <- length(lambda)
    llVec <- rep(NA, lenlam)
    for (i in 1:lenlam)
    {
#        lv <- lambda[i]
#        assign("bcFct", bfFct(lambda[i]), envir = .GlobalEnv)     
        assign("bcFct1", bfFct(lambda[i]), envir = .GlobalEnv)     
        assign("bcFct2", bfFct2(lambda[i]), envir = .GlobalEnv)     
        newFormula <- bcFct1(.) ~ bcFct2(.)        

        nlsTemp <- try(update(object, formula. = newFormula, start = startVec, trace = FALSE), silent = TRUE)
#        if (!inherits(nlsTemp, "try-error")) {llVec[i] <- logLik(nlsTemp)}
        if (!inherits(nlsTemp, "try-error")) {llVec[i] <- llFct(nlsTemp, lambda[i])}
    }

    lv <- lambda[which.max(llVec)]
    llv <- max(llVec, na.rm = TRUE)
    ci <- boxcoxCI(lambda, llVec, level)
    
    if (plotit)  # based on boxcox.default
    {
        plot(lambda, llVec, type ="l", xlab = xlab, ylab = ylab, ...)
        
        plims <- par("usr")
        y0 <- plims[3]
        lim <- llv-qchisq(level, 1)/2
        
        segments(lv, llv, lv, y0, lty=3)
        segments(ci[1], lim, ci[1], y0, lty = 3)  # lower limit
        segments(ci[2], lim, ci[2], y0, lty = 3)  # upper limit
        
        scal <- (1/10 * (plims[4] - y0))/par("pin")[2] 
        scx <- (1/10 * (plims[2] - plims[1]))/par("pin")[1] 
        text(lambda[1] + scx, lim + scal, " 95%") 
        abline(h = lim, lty = 3)

    }
#    assign("bcFct", bfFct(lv), envir = .GlobalEnv) 
    assign("bcFct1", bfFct(lv), envir = .GlobalEnv) 
    assign("bcFct2", bfFct2(lv), envir = .GlobalEnv) 
    retFit <- update(object, formula. = newFormula, start = startVec, trace = FALSE)  
    # last time 'bcFct1' and 'bcFct2' are used
#    rm(bcFct1, bcFct2, envir = .GlobalEnv)  # to ensure that the predict method can find bcFct2()
    retFit$lambda <- list(lambda = lv, ci = ci)
    invisible(retFit)
}


"boxcoxCI" <- 
function(x, y, level = 0.95)
{
    ## R lines taken from boxcox.default in the package MASS and then slightly modified
    xl <- x
    loglik <- y
    
    llnotna <- !is.na(loglik)
    xl <- xl[llnotna]
    loglik <- loglik[llnotna]
    
    m <- length(loglik)

    mx <- (1:m)[loglik == max(loglik)][1]
    Lmax <- loglik[mx]
    lim <- Lmax - qchisq(level, 1)/2

    ind <- range((1:m)[loglik > lim])

    xx <- rep(NA, 2)
    if(loglik[1] < lim) 
    {
        i <- ind[1]
        xx[1] <- xl[i - 1] + ((lim - loglik[i - 1]) *
                          (xl[i] - xl[i - 1]))/(loglik[i] - loglik[i - 1])

    }
    if(loglik[m] < lim) 
    {
        i <- ind[2] + 1
        xx[2] <- xl[i - 1] + ((lim - loglik[i - 1]) *
                          (xl[i] - xl[i - 1]))/(loglik[i] - loglik[i - 1])
    }
    return(xx)
}

"bcSummary" <- function(object)
{
    cat("\n")
    cat("Estimated lambda:", object$lambda$lambda, "\n")
 
    bcci <- format(object$lambda$ci, digits = 2)
    ciStr <- paste("[", bcci[1], ",", bcci[2], "]", sep="")
    cat("Confidence interval for lambda:", ciStr, "\n\n")    
}

