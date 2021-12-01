#Asymptotic regression model ##############################
asymReg.fun <- function(predictor, init, m, plateau) {
                      x <- predictor
                      plateau - (plateau - init) * exp (- m * x)
}

asymReg.Init <- function(mCall, LHS, data, ...) {
          xy <- sortedXyData(mCall[["predictor"]], LHS, data)
          x <-  xy[, "x"]; y <- xy[, "y"]
          plateau <- NLSstRtAsymptote(xy)
          dataR <- xy[1:2,]
          init <- coef(lm(y ~ x, data = dataR))[1]
      
          ## Linear regression on pseudo y values
          pseudoY <- log( ( y - plateau )/(init - plateau ) )
          coefs <- coef( lm(pseudoY ~ x - 1) )
          m <- - coefs[1]
          value <- c(init, m, plateau)
          names(value) <- mCall[c("init", "m", "plateau")]
          value
}

NLS.asymReg <- selfStart(asymReg.fun, asymReg.Init, parameters=c("init", "m", "plateau"))

DRC.asymReg <- function(fixed = c(NA, NA, NA), names = c("init", "m", "plateau")){
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

        init <- parmMat[, 1]; m <- parmMat[, 2]; plateau <- parmMat[, 3]; 
        plateau - (plateau - init) * exp (- m * x)
    }

    ## Defining self starter function
    ssfct <- function(dataf)
    {
        x <- dataf[, 1]
        y <- dataf[, 2]

        dataS <- sortedXyData(x, y)
        plateau <- NLSstRtAsymptote(dataS)
        dataR <- dataS[1:2,]
        init <- coef(lm(y ~ x, data = dataR))[1]
      
        ## Linear regression on pseudo y values
        pseudoY <- log( ( y - plateau )/(init - plateau ) )
        coefs <- coef( lm(pseudoY ~ x - 1) )
        m <- - coefs[1]
        return(c(init, m, plateau)[notFixed])
    }

    ## Defining names
    pnames <- names[notFixed]

    ## Defining derivatives

    ## Defining the ED function

    ## Defining the inverse function

    ## Defining descriptive text
    text <- "Asymptotic Regression Model"

    ## Returning the function with self starter and names
    returnList <- list(fct = fct, ssfct = ssfct, names = pnames, text = text, noParm = sum(is.na(fixed)))

    class(returnList) <- "drcMean"
    invisible(returnList)
}

"DRC.SSasymp" <- function(fixed = c(NA, NA, NA), names = c("Asym", "R0", "lrc")) {
    ## Checking arguments
    numParm <- 3
    
    # if (!is.character(names) | !(length(names) == numParm)) {stop("Not correct 'names' argument")}
    # if (!(length(fixed) == numParm)) {stop("Not correct 'fixed' argument")}

    ## Fixing parameters (using argument 'fixed')
    notFixed <- is.na(fixed)
    parmVec <- rep(0, numParm)
    parmVec[!notFixed] <- fixed[!notFixed]

    ## Defining the non-linear function
    fct <- function(x, parm)
    {
        parmMat <- matrix(parmVec, nrow(parm), numParm, byrow = TRUE)
        parmMat[, notFixed] <- parm

        Asym <- parmMat[, 1]; R0 <- parmMat[, 2]; lrc <- parmMat[, 3]
        Asym + (R0 - Asym) * exp (- exp(lrc) * x)
    }

    ## Defining self starter function
    ssfct <- function(dataf)
    {
      
      x <- dataf[, 1]
      y <- dataf[, 2]
      dataS <- sortedXyData(x, y)
      Asym <- NLSstRtAsymptote(dataS)
      dataR <- dataS[1:2,]
      R0 <- coef(lm(y ~ x, data = dataR))[1]
      
        ## Linear regression on pseudo y values
        pseudoY <- log( ( y - Asym )/(R0 - Asym ) )
        coefs <- coef( lm(pseudoY ~ x - 1) )
        b <- - coefs[1]
        lrc <- log(b)
        return(c(Asym, R0, lrc)[notFixed])
    }

    ## Defining names
    pnames <- names[notFixed]

    ## Defining derivatives

    ## Defining the ED function

    ## Defining the inverse function

    ## Defining descriptive text
    text <- "Asymptotic regression model"

    ## Returning the function with self starter and names
    returnList <- list(fct = fct, ssfct = ssfct, names = pnames, text = text, noParm = sum(is.na(fixed)))

    class(returnList) <- "drcMean"
    invisible(returnList)
}

# #Growth
# set.seed(1234)
# X <- c(1, 3, 5, 7, 9, 11, 13, 20)
# plateau <- 20; init <- 5; m <- 0.3
# Ye <- asymReg.fun(X, init, m, plateau)
# epsilon <- rnorm(8, 0, 0.5)
# Y <- Ye + epsilon
# model <- drm(Y ~ X, fct = DRC.asymReg())
# plot(model, log="", main = "Asymptotic regression")
# model2 <- nls(Y ~ NLS.asymReg(X, init, m, plateau))
# model3 <- nls(Y ~ SSasymp(X, Asymp, R0, lrc))
# model4 <- drm(Y ~ X, fct = DRC.SSasymp())
# plot(model, log="", main = "Asymptotic regression")
# summary(model)
# summary(model2)
# summary(model3)
# exp(-1.0980)
# summary(model4)
# 
# set.seed(1234)
# X <- c(1, 3, 5, 7, 9, 11, 13, 20)
# plateau <- 5; init <- 20; m <- -1.08
# Ye <- SSasymp(X, plateau, init, m)
# epsilon <- rnorm(8, 0, 0.5)
# Y <- Ye + epsilon
# model <- drm(Y ~ X, fct = DRC.SSasymp())
# plot(model, log="", main = "Asymptotic regression")
# model2 <- nls(Y ~ SSasymp(X, Asymp, R0, lrc))
# summary(model)
# summary(model2)
