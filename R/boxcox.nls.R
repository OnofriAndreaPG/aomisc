"boxcox.nls" <- function(object, lambda = seq(-2, 2, 1/10), plotit = TRUE, start,
eps = 1/50, bcAdd = 0, level = 0.95, xlab = expression(lambda), ylab = "log-likelihood", ...)
{
  
  # Defining the transformation function
  bfFct <- function(x, lv, eps, bcAdd){
    if (abs(lv) < eps) {return(log(x+bcAdd))} else {return(((x+bcAdd)^lv-1)/lv)}
  }
  
  # Defining the new fitting function
  newFct <- function(x, lv, eps, bcAdd){
    if (abs(lv) < eps) {
      return(paste("log(", x, ")", sep = ""))
    } else {
      .tmp <- paste("(", x, "+", bcAdd, ")", sep = "")
      .tmp <- paste(.tmp, "^(", lv, ")", sep = "")
      .tmp <- paste("(", .tmp, "- 1)", sep = "")
      .tmp <- paste(.tmp, "/(", lv, ")", sep = "")
      return(.tmp)
    }
  }
  
  # Retrieving objects from model call
  df <- eval(object$data)
  mm <- object$m
  cc <- object$call
  pnms <- names(mm$getPars())
  form <- cc$formula
  rhsnms <- all.vars(form[[3]])
  vnms <- rhsnms[!(rhsnms %in% pnms)]
  namList <- list(x = as.name(vnms), y = form[[2]])
  x <- df[,as.character(namList$x)]
  y <- df[,as.character(namList$y)]
  oldFormula <- as.character(formula(object))[[3]]
  
  # Re-fitting the model with TBS approach (recursive)
  lenlam <- length(lambda)
  llVec <- rep(NA, lenlam)
  rss <- rep(NA, lenlam)
  k <- rep(NA, lenlam)
  for (i in 1:lenlam) {
    df$newRes <- bfFct(y, lambda[i], eps, bcAdd)
    newFormula <- newFct(oldFormula, lambda[i], eps, bcAdd)
    nlsTemp <- try(nls(formula = newRes ~ eval(parse(text = newFormula), df), 
                  data = df, 
                  start = coef(object)), silent = TRUE)
    if (!inherits(nlsTemp, "try-error")) {
      # Calculate log-likelihood
      sumObj <- summary(nlsTemp)
      rss[i] <- deviance(nlsTemp)
      N <- sum(sumObj$df)
      Ji <- y^(lambda[i] - 1)
      llVec[i] <- -N * log(sqrt(sum(residuals(nlsTemp)^2)/N))-N/2+sum(log(Ji))
    } else {
      llVec[i] <- NA
      rss[i] <- NA
    }
  }
  
  lv <- lambda[which.max(llVec)]
  llv <- max(llVec, na.rm = TRUE)
  rssF <- rss[which.max(llVec)]
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
  
  # Refit model based on best lambda
  df$newRes <- bfFct(y, lv, eps, bcAdd)
  newFormula <- newFct(oldFormula, lv, eps, bcAdd)
  retFit <- try(nls(formula = newRes ~ eval(parse(text = newFormula), df), 
                  data = df, 
                  start = coef(object)), silent = TRUE)
  retFit$lambda <- list(lambda = lv, ci = ci, loglik = llv)
  class(retFit) <- c("nlsbc", "nls")
  invisible(retFit)
  # return(returnList)
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

"summary.nlsbc" <- function(object, ...)
{
    cat("\n")
    cat("Estimated lambda:", object$lambda$lambda, "\n")
 
    bcci <- format(object$lambda$ci, digits = 2)
    ciStr <- paste("[", bcci[1], ",", bcci[2], "]", sep="")
    cat("Confidence interval for lambda:", ciStr, "\n\n")
    class(object) <- "nls"
    summary(object)
}

# "boxcox.nlsOLD" <- function(object, lambda = seq(-2, 2, 1/10), plotit = TRUE, start,
# eps = 1/50, bcAdd = 0, level = 0.95, xlab = expression(lambda), ylab = "log-likelihood", ...)
# {
#     ## Defining the Box-Cox modified power transformations
#     bfFct <- function(lv)
#     {
#         function(x) {if (abs(lv) < eps) {return(log(x+bcAdd))} else {return(((x+bcAdd)^lv-1)/lv)}}    
#     }
# 
#     evalForm <- eval(formula(object)[[3]][[1]])
#     bfFct2 <- function(lv)
#     {
#         bcFct2 <- function(x)
#         {
#             ## Transforming the mean
#             bcFct <- bfFct(lv)
#             bcFctVal <- bcFct(x)
# 
# #            ## Assigning sttributes
# #            attrList <- attributes(evalForm)
# #            lenAL <- length(attrList)
# #            namAL <- names(attrList)
# #     
# #            for (i in 1:lenAL)
# #            {
# ##                attr(bcFctVal, namAL[i]) <- attrList[[i]]
# #            }
# #            
#             ## Adjusting gradient
#             grad1 <- attr(bcFctVal, "gradient")
#             grad2 <- grad1*(x^(lv-1))
#             attr(bcFctVal, "gradient") <- grad2 
#             
#             bcFctVal
#         }
#         bcFct2
#      }
#      
#      if (!inherits(evalForm , "selfStart"))
#      {
#          bfFct2 <- bfFct
#      }    
#          
# #    bcFct <- function(x) {if (abs(lv) < eps) {return(log(x+bcAdd))} else {return(((x+bcAdd)^lv-1)/lv)}}
# #    assign("bcFct", bcFct, envir = .GlobalEnv)
# #    newFormula <- bcFct(.) ~ bcFct(.)
#     if (missing(start)) 
#     {
#         startVec <- coef(object)
#     } else {
#         startVec <- start
#     }
#     
#     ## Defining the log-likelihood function
#     llFct <- function(object, lv)
#     {
#         # sumObj <- summary(object)
#         N <- length(residuals(object)) # sum(sumObj$df)
#         Ji <- (eval(object$data)[, as.character(formula(object)[[2]])[2]])^(lv-1)
#         -N*log(sqrt(sum(residuals(object)^2)/N))-N/2+sum(log(Ji))
#     }
#     
#     lenlam <- length(lambda)
#     llVec <- rep(NA, lenlam)
#     rss <- rep(NA, lenlam)
#     k <- rep(NA, lenlam)
#     # coefs <- data.frame(uno, due)
#     for (i in 1:lenlam)
#     {
# #        lv <- lambda[i]
# #        assign("bcFct", bfFct(lambda[i]), envir = .GlobalEnv)     
#         assign("bcFct1", bfFct(lambda[i]), envir = .GlobalEnv)     
#         assign("bcFct2", bfFct2(lambda[i]), envir = .GlobalEnv)     
#         newFormula <- bcFct1(.) ~ bcFct2(.)        
#         
#         nlsTemp <- try(update(object, formula. = newFormula, start = startVec, trace = FALSE), silent = TRUE)
# #        if (!inherits(nlsTemp, "try-error")) {llVec[i] <- logLik(nlsTemp)}
#         if (!inherits(nlsTemp, "try-error")) {
#           llVec[i] <- llFct(nlsTemp, lambda[i])
#           rss[i] <- deviance(nlsTemp)
#           k[i] <- coef(nlsTemp)[2]
#         }
#         # print(coef(nlsTemp))
#     }
# 
#     lv <- lambda[which.max(llVec)]
#     llv <- max(llVec, na.rm = TRUE)
#     ci <- boxcoxCI(lambda, llVec, level)
#     
#     if (plotit)  # based on boxcox.default
#     {
#         plot(lambda, llVec, type ="l", xlab = xlab, ylab = ylab, ...)
#         
#         plims <- par("usr")
#         y0 <- plims[3]
#         lim <- llv-qchisq(level, 1)/2
#         
#         segments(lv, llv, lv, y0, lty=3)
#         segments(ci[1], lim, ci[1], y0, lty = 3)  # lower limit
#         segments(ci[2], lim, ci[2], y0, lty = 3)  # upper limit
#         
#         scal <- (1/10 * (plims[4] - y0))/par("pin")[2] 
#         scx <- (1/10 * (plims[2] - plims[1]))/par("pin")[1] 
#         text(lambda[1] + scx, lim + scal, " 95%") 
#         abline(h = lim, lty = 3)
# 
#     }
# #    assign("bcFct", bfFct(lv), envir = .GlobalEnv) 
#     assign("bcFct1", bfFct(lv), envir = .GlobalEnv) 
#     assign("bcFct2", bfFct2(lv), envir = .GlobalEnv) 
#     retFit <- update(object, formula. = newFormula, start = startVec, trace = FALSE)  
#     # last time 'bcFct1' and 'bcFct2' are used
# #    rm(bcFct1, bcFct2, envir = .GlobalEnv)  # to ensure that the predict method can find bcFct2()
#     retFit$lambda <- list(lambda = lv, ci = ci, llVec = llVec, rss = rss, k = k)
#     invisible(retFit)
# }


# "boxcoxCI" <- 
# function(x, y, level = 0.95)
# {
#     ## R lines taken from boxcox.default in the package MASS and then slightly modified
#     xl <- x
#     loglik <- y
#     
#     llnotna <- !is.na(loglik)
#     xl <- xl[llnotna]
#     loglik <- loglik[llnotna]
#     
#     m <- length(loglik)
# 
#     mx <- (1:m)[loglik == max(loglik)][1]
#     Lmax <- loglik[mx]
#     lim <- Lmax - qchisq(level, 1)/2
# 
#     ind <- range((1:m)[loglik > lim])
# 
#     xx <- rep(NA, 2)
#     if(loglik[1] < lim) 
#     {
#         i <- ind[1]
#         xx[1] <- xl[i - 1] + ((lim - loglik[i - 1]) *
#                           (xl[i] - xl[i - 1]))/(loglik[i] - loglik[i - 1])
# 
#     }
#     if(loglik[m] < lim) 
#     {
#         i <- ind[2] + 1
#         xx[2] <- xl[i - 1] + ((lim - loglik[i - 1]) *
#                           (xl[i] - xl[i - 1]))/(loglik[i] - loglik[i - 1])
#     }
#     return(xx)
# }
# 
# "bcSummary" <- function(object)
# {
#     cat("\n")
#     cat("Estimated lambda:", object$lambda$lambda, "\n")
#  
#     bcci <- format(object$lambda$ci, digits = 2)
#     ciStr <- paste("[", bcci[1], ",", bcci[2], "]", sep="")
#     cat("Confidence interval for lambda:", ciStr, "\n\n")    
# }

