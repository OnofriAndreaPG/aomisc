getTau <- function(formula, vi, data, method="HE"){
  # This function calculates an estimate for tau2 in meta-analyses
  # Updated on 10/06/2022
  MLfun <- function(X, Y, vi, tau){
    totVar <- tau + vi
    w <- 1/totVar
    modlmW <- lm.wfit(X, Y, w)
    - ( -length(Y)/2*log(2*pi)-1/2*sum(log(totVar))-
          1/2*sum(modlmW$residuals^2/totVar ) )
  }
  REMLfun <- function(X, Y, vi, tau){
    totVar <- tau + vi
    w <- 1/totVar
    modlmW <- lm.wfit(X, Y, w)

    - ( -length(Y)/2*log(2*pi)-1/2*sum(log(totVar))-
          1/2*sum(modlmW$residuals^2/totVar ) -
          0.5*log(sum(1/totVar)))
  }
  
  callDetail <- match.call()

  ## Handling the 'formula', 'vi' and 'data' arguments ##########
  anName <- deparse(substitute(vi))  # storing name for later use
  if (length(anName) > 1) {anName <- anName[1]}  # to circumvent the behaviour of 'substitute' in do.call("multdrc", ...)
  if (nchar(anName) < 1) {stop("no standad errors are provided")}  # in case only one curve is analysed

  mf <- match.call(expand.dots = FALSE)
  nmf <- names(mf)
  mnmf <- match(c("formula", "vi", "data"), nmf, 0)

  mf[[1]] <- as.name("model.frame")

  mf <- eval(mf[c(1, mnmf)], parent.frame())  #, globalenv())
  mt <- attr(mf, "terms")
  varNames <- names(mf)[c(2, 1)]  # Rigido x1 + y
  varNames0 <- names(mf) # tutte le variabili
  viVals <- model.extract(mf, "vi")


  # only used once, but mf is overwritten later on
  dose <- model.matrix(mt, mf)[,-c(1)]  # with no intercept
  xDim <- ncol(as.matrix(dose))
  Y <- model.response(mf, "numeric")
  
  X <- model.matrix(~dose)
  naiveMod <- lm.fit(X, Y)
  
  
  # MSE <- sum(naiveMod$residuals^2)/naiveMod$df
  # tau <- MSE - mean(vi) Va bene solo per i casi di ANOVA models
  RSS <- sum(naiveMod$residuals^2)
  
  k <- length(Y)
  p <- ncol(X)
  H <- X %*% solve(t(X)%*%X) %*% t(X) 
  P     <- diag(k) -  H # p is the central matrix A
  V     <- diag(viVals, nrow=k, ncol=k)
  PV    <- P %*% V # note: this is not symmetric
  trPV  <- sum(diag(PV)) # since PV needs to be computed anyway, can use 
  tau  <- (RSS - trPV) / (k-p)

  if(method!="HE"){
    if(method=="ML"){
      likfun <- optim(par = tau, lower = 0,
                      fn = MLfun, hessian = T,
                      method = "L-BFGS-B",
                      Y = Y, vi = viVals, X=X)
      tau <- likfun$par
      # print(sqrt(diag(solve(likfun$hessian))))
      }
    else { if(method=="REML") {
      likfun <- optim(par = tau, lower = 0,
                      fn = REMLfun,
                      method = "L-BFGS-B",
                      Y = Y, vi = viVals, X = X)
      tau <- likfun$par }
    }}
  tau}

myMeta <- function(Y, vi, mods, data, method="HE"){
  tau <- getTau(Y, vi, mods, data, method)
  totVar <- tau + vi
  w <- 1/totVar
  W <- diag(w)
  X <- model.matrix(mods, data = data)
  modlmW <- lm.wfit(X, Y, w)
  #summary(modlmW)
  #deviance(modlmW)
  #residuals(modlmW)
  #wls <- sum(w*(modlmW$residuals^2)) #This is minimised
  #(resid_var <- sum(wls/modlmW$df.residual) )
  #sqrt(resid_var) #Residual standard deviation
  #sqrt(1/sum(1/totVar)) #SE
  n <- length(Y)
  terms <- modlmW$coef
  vcovMA <- solve(t(X) %*% W %*% X) %*%
    (t(X) %*% W %*% diag(totVar, n) %*%
       t(W) %*% X) %*% solve(t(X) %*% W %*% X)
  es <- sqrt(diag(vcovMA))
  list("tau" = tau, "coeff" = data.frame("Estimate"=terms, "SE"= es), "vcov" = vcovMA)
}

# #Prova
# prova <- function(formula, fct, data){
#   fr <- model.frame(formula, data)
#   X <- model.matrix(fr, data)
#   Y <- model.response(fr, "numeric")
#   mod <-  drm(formula, fct=fct, type="event", data=data)
#   summary(mod)$coef
# }
# #Prova
# # verbascum <- read.csv("ThreeGenotypes.csv", header=T)
# # datino <- subset(verbascum, Dish==265)
# # prova(nSeeds ~ timeBef + timeAf, fct=LL.3(), data=datino)
# #
# # formula <- ED50 ~ Temp
# # data <- step2Data
# # fct = linear()

# getTauDrc <- function(formula, SE, fct, data, method="HE"){
#   MLfun <- function(formula, SE, fct, data, parms){
#     fr <- model.frame(formula, data)
#     X <- model.matrix(fr, data)[,2]
#     Y <- model.response(fr, "numeric")
#     vi <- SE^2
#     totVar <- parms[1] + vi
#     w <- 1/totVar
#     parmMat <- matrix(parms[2:length(parms)], 1, length(parms[2:length(parms)]))
#     expt <- fct$fct(X, parmMat)
#     res <- Y - expt
#     - ( -length(Y)/2*log(2*pi)-1/2*sum(log(totVar))-
#           1/2*sum(res^2/totVar ) )
#   }
#   REMLfun <- function(formula, SE, fct, data, parms){
#     fr <- model.frame(formula, data)
#     X <- model.matrix(fr, data)[,2]
#     Y <- model.response(fr, "numeric")
#     vi <- SE^2
#     totVar <- parms[1] + vi
#     w <- 1/totVar
#     parmMat <- matrix(parms[2:length(parms)], 1, length(parms[2:length(parms)]))
#     expt <- fct$fct(X, parmMat)
#     res <- Y - expt
#     - ( -length(Y)/2*log(2*pi)-1/2*sum(log(totVar))-
#           1/2*sum(res^2/totVar ) -
#           0.5*log(sum(1/totVar)))
#   }
#   #data=step2Data
#   #formula <- ED50 ~ Temp
#   #fct <- linear()
#   SE <- SE
#   vi <- SE^2
#   naiveMod <- drm(formula, data=data, fct=fct)
#   MSE <- sum(residuals(naiveMod)^2)/naiveMod$df
#   tau <- MSE - mean(vi)
#   if(method!="HE"){
#     parms <- c(tau, naiveMod$coef)
#     if(method=="ML"){
#       likfun <- optim(par = parms, lower = c(0, rep(NA, length(parms))),
#                       fn = MLfun,
#                       method = "L-BFGS-B", hessian=T,
#                       formula = formula, SE = SE, fct = fct, data=data)
#       tau <- likfun }
#     else { if(method=="REML") {
#       naiveMod$coefs
#       likfun <- optim(par = parms, lower = c(0, rep(NA, length(parms))),
#                       fn = REMLfun,
#                       method = "L-BFGS-B", hessian=T,
#                       formula = formula, SE = SE, fct = fct, data=data)
#       tau <- likfun }
#     }}
#   tau}
#
# metaDrm <- function(formula, SE, fct, data, method="HE"){
#   tau <- getTauDrc(formula, SE, fct, data, method)
#   fr <- model.frame(formula, data)
#   X <- model.matrix(fr, data)
#   Y <- model.response(fr, "numeric")
#   vi <- data$SE^2
#   totVar <- tau + vi
#   w <- 1/totVar
#   W <- diag(w)
#   #modlmW <- lm.wfit(X, Y, w)
#   modlmW <- drm(formula, weights=sqrt(w),
#                 data=data, fct=fct)
#   #summary(modlmW)
#   #deviance(modlmW)
#   #residuals(modlmW)
#   #wls <- sum(w*(modlmW$residuals^2)) #This is minimised
#   #(resid_var <- sum(wls/modlmW$df.residual) )
#   #sqrt(resid_var) #Residual standard deviation
#   #sqrt(1/sum(1/totVar)) #SE
#   n <- length(Y)
#   terms <- modlmW$coef
#   vcovMA <- solve(t(X) %*% W %*% X) %*%
#     (t(X) %*% W %*% diag(totVar, n) %*%
#        t(W) %*% X) %*% solve(t(X) %*% W %*% X)
#   es <- sqrt(diag(vcovMA))
#   list("tau" = tau, "coeff" = data.frame("Estimate"=terms, "SE"= es), "vcov" = vcovMA)
# }
#
# getTauDrc <- function(formula, SE, fct, data, method="HE"){
#   MLfun <- function(formula, SE, fct, data, parms){
#     fr <- model.frame(formula, data)
#     X <- model.matrix(fr, data)[,2]
#     Y <- model.response(fr, "numeric")
#     vi <- SE^2
#     totVar <- parms[1] + vi
#     w <- 1/totVar
#     parmMat <- matrix(parms[2:length(parms)], 1, length(parms[2:length(parms)]))
#     expt <- fct$fct(X, parmMat)
#     res <- Y - expt
#     - ( -length(Y)/2*log(2*pi)-1/2*sum(log(totVar))-
#           1/2*sum(res^2/totVar ) )
#   }
#   REMLfun <- function(formula, SE, fct, data, parms){
#     fr <- model.frame(formula, data)
#     X <- model.matrix(fr, data)[,2]
#     Y <- model.response(fr, "numeric")
#     vi <- SE^2
#     totVar <- parms[1] + vi
#     w <- 1/totVar
#     parmMat <- matrix(parms[2:length(parms)], 1, length(parms[2:length(parms)]))
#     expt <- fct$fct(X, parmMat)
#     res <- Y - expt
#     - ( -length(Y)/2*log(2*pi)-1/2*sum(log(totVar))-
#           1/2*sum(res^2/totVar ) -
#           0.5*log(sum(1/totVar)))
#   }
#   #data=step2Data
#   #formula <- ED50 ~ Temp
#   #fct <- linear()
#   SE <- SE
#   vi <- SE^2
#   naiveMod <- drm(formula, data=data, fct=fct)
#   MSE <- sum(residuals(naiveMod)^2)/naiveMod$df
#   tau <- MSE - mean(vi)
#   if(method!="HE"){
#     parms <- c(tau, naiveMod$coef)
#     if(method=="ML"){
#       likfun <- optim(par = parms, lower = c(0, rep(NA, length(parms))),
#                       fn = MLfun,
#                       method = "L-BFGS-B", hessian=T,
#                       formula = formula, SE = SE, fct = fct, data=data)
#       tau <- likfun }
#     else { if(method=="REML") {
#       naiveMod$coefs
#       likfun <- optim(par = parms, lower = c(0, rep(NA, length(parms))),
#                       fn = REMLfun,
#                       method = "L-BFGS-B", hessian=T,
#                       formula = formula, SE = SE, fct = fct, data=data)
#       tau <- likfun }
#     }}
#   tau}


# Funzione per il calcolo solo con MLE###############################

# #This is the original function in metafor (not very transparent!)
# getTau <- function(Y, vi, X, method="HE"){
#   change <- 10
#   maxiter <- 100
#   threshold <- 10E-05
#   iter <- 0
#   k <- length(Y)
#   tau2.min <- 0
#   stepadj = 1
#   p <- length(X[1,])
#   #Method HE
#   stXX <- solve(t(X) %*% X)
#   P     <- diag(k) - X %*% tcrossprod(stXX,X)
#   resid <- crossprod(Y, P)
#   RSS   <- resid %*% Y #Calcolo matriciale della RSS di un modello naive
#   #RSS <- resid %*% resid
#   V     <- diag(vi, nrow=k, ncol=k)
#   PV    <- P %*% V ### careful: is not symmetric => restituisce il divisore di 7/8
#   trPV  <- sum(diag(PV)) ### since PV needs to be computed anyway, can use .tr()
#   tau2  <- (RSS - trPV) / (k-p) #Differenze tra devianza del modello naive e devianza media dei SES
#   tau2 <- as.vector(tau2)
#
#   if(method!="HE"){
#     while (change > threshold) {
#       iter <- iter + 1
#       old2 <- as.vector(tau2)
#       tau2 <- as.vector(tau2)
#       wi   <- 1/(vi + tau2)
#       if (any(tau2 + vi < 0))
#         stop(mstyle$stop("Some marginal variances are negative."))
#       if (any(is.infinite(wi)))
#         stop(mstyle$stop("Division by zero when computing the inverse variance weights."))
#       W <- diag(wi, nrow=k, ncol=k)
#       #stXWX <- .invcalc(X=X, W=W, k=k) #Funzione sottostante
#       stXWX <- solve(t(X) %*% W %*% X)
#       P   <- W - W %*% X %*% stXWX %*% crossprod(X,W)
#
#       if (method == "ML") {
#         PP  <- P %*% P
#         adj <- (crossprod(Y,PP) %*% Y - sum(wi)) / sum(wi^2)
#       }
#       if (method == "REML") {
#         PP  <- P %*% P
#         adj <- (crossprod(Y, PP) %*% Y - sum(diag(P))) / sum(diag(PP))
#       }
#       if (method == "EB") {
#         adj <- (crossprod(Y,P) %*% Y * k/(k-p) - k) / sum(wi)
#       }
#
#       adj <- adj * stepadj ### apply (user-defined) step adjustment
#
#       while (tau2 + adj < tau2.min) ### use step-halving if necessary
#         adj <- adj / 2
#
#       #tau2   <- ifelse(tau2.fix, tau2.val, tau2 + adj)
#       tau2   <- tau2 + adj
#       change <- abs(old2 - tau2)
#
#       if (iter > maxiter) {
#         conv <- 0
#         break
#       }
#
#     }
#   }
#   as.vector(tau2)}
#
