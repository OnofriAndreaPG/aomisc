gnlht <- function(object, ...){
  # Edited on 12/1/21 #
   UseMethod("gnlht")
}
gnlht.old <- function(obj, func,  const = NULL, vcov. = vcov, parameterNames = names(coef(obj))){
   
   temp <- lapply(func, function(x) as.character(as.expression(x[[length(x)]])))
   func <- data.frame(form = unlist(temp))
   if(any(class(obj) == "nlme") | any(class(obj) == "lme")) { 
     coefs <- fixef(obj)
   } else { coefs <- coef(obj) }
   names(coefs) <- parameterNames
   
   if (is.function(vcov.)){
            vcMat <- vcov.(obj)
        } else {vcMat <- vcov.}
   
 
   if(is.null(const) == T){
      inputs <- func
      colnames(inputs) <- c("Form") 
      val <- data.frame()
      for(i in 1:length(inputs[,1])){
      res <- car::deltaMethod(object=coefs, g=as.character(inputs$Form[i]), 
                      vcov.=vcMat, level = 0.95 )
      val <- rbind(val, res)
     }
     row.names(val) <- 1:length(inputs[,1])
     val <- cbind(inputs, val)
     # lisRes <- apply(func, 1, 
     #               function(x) car::deltaMethod(object=coefs, g=x, 
     #                    vcov.=vcMat, level = 0.95))
     # val <- plyr::ldply(lisRes)
     # val <- cbind(func, val)
   }else{
      # Una sola costante è possibile, anche se può assumere più valori
     inputs <- expand.grid(func$form, const[,1])
     inputs
     colnames(inputs) <- c("Form", colnames(const)[1]) 
     val <- data.frame()
     for(i in 1:length(inputs[,1])){
     valCost <- inputs[i, 2]
     names(valCost) <- colnames(inputs)[2]
     res <- car::deltaMethod(object=coefs, g=as.character(inputs$Form[i]), 
                      vcov.=vcMat, level = 0.95, 
                      constants = valCost)
     val <- rbind(val, res)
     }
     row.names(val) <- 1:length(inputs[,1])
     val <- cbind(inputs, val)
     # val
     # Corrections made on 4/11/2020 - Andrea
     # lisRes <- apply(const, 1, function(y) apply(func, 1, 
     #               function(x) car::deltaMethod(object=coefs, g=x, 
     #                    vcov.=vcMat, level = 0.95, constants = y )))
     # 
     # val <- plyr::ldply(do.call(rbind, lisRes))
     # funList <- func %>% slice(rep(1:n(), each = length(const[,1])))
     # func
     # constList <- const %>% slice(rep(1:n(), length(func[,1])))
     # val <- cbind(funList, constList, val)
   }
  #val
  lenVal <- length(val[1,]) 
  retDF <- val[, c(-(lenVal-1), -lenVal)]
  retDF$"t-value" <- abs(retDF$Estimate/retDF$SE)
  if(class(obj)[1] == "nls") resDF <- summary(obj)$df[2]
  else if(class(obj)[1] == "drc") resDF <- summary(obj)$df
  else if(class(obj)[1] == "lm") resDF <- obj$df.residual
  else resDF <- Inf
  retDF$"p-value" <- 2 * pt(retDF$"t-value", resDF, lower.tail = F)
  if(resDF == Inf) colnames(retDF)[length(colnames(retDF)) - 1] <- "Z-value"
  return(retDF)
}
gnlht.nlme <- function(object, func,  const = NULL, vcov. = vcov, 
                     parameterNames = NULL, dfr = NULL){
   coefs <- fixef(object)
   dfr <- summary(modNlin)$df[2]   
   if (is.function(vcov.)){
            vcMat <- vcov.(object)
        } else {vcMat <- vcov(object)}
   
  retList <- gnlht.default(coefs = coefs, func = func,  const = const,
                           vcov. = vcMat, parameterNames = names(coef(object)),
                           dfr = dfr)
  return(retList)
}   

gnlht.lme <- function(object, func,  const = NULL, vcov. = vcov, 
                     parameterNames = NULL, dfr = NULL){
   coefs <- fixef(object)
   dfr <- summary(modNlin)$df[2]   
   if (is.function(vcov.)){
            vcMat <- vcov.(object)
        } else {vcMat <- vcov(object)}
   
  retList <- gnlht.default(coefs = coefs, func = func,  const = const,
                           vcov. = vcMat, parameterNames = names(coef(object)),
                           dfr = dfr)
  return(retList)
}   

gnlht.lm <- function(object, func,  const = NULL, vcov. = vcov, 
                     parameterNames = NULL, dfr = NULL){
   coefs <- coef(object)
   dfr <- summary(modNlin)$df[2]   
   if (is.function(vcov.)){
            vcMat <- vcov.(object)
        } else {vcMat <- vcov(object)}
   
  retList <- gnlht.default(coefs = coefs, func = func,  const = const,
                           vcov. = vcMat, parameterNames = names(coef(object)),
                           dfr = dfr)
  return(retList)
}   

gnlht.nls <- function(object, func,  const = NULL, vcov. = vcov, parameterNames = NULL, dfr = NULL){
   coefs <- coef(object)
   dfr <- summary(modNlin)$df[2]   
   if (is.function(vcov.)){
            vcMat <- vcov.(object)
        } else {vcMat <- vcov(object)}
   
  retList <- gnlht.default(coefs = coefs, func = func,  const = const,
                           vcov. = vcMat, parameterNames = names(coef(object)),
                           dfr = dfr)
  return(retList)
}   
   
gnlht.numeric <- function(object, func,  const = NULL, vcov. = vcov, parameterNames = NULL, dfr = NULL){

  retList <- gnlht.default(coefs = object, func = func,  const = const,
                           vcov. = vcov., parameterNames = parameterNames,
                           dfr = dfr)
  return(retList)
}   

gnlht.default <- function(coefs, func,  const = NULL, vcov. = vcov, parameterNames = NULL, dfr = NULL){
   
   temp <- lapply(func, function(x) as.character(as.expression(x[[length(x)]])))
   func <- data.frame(form = unlist(temp))
   names(coefs) <- parameterNames
   vcMat <- vcov.
   
   if(is.null(const) == T){
      inputs <- func
      colnames(inputs) <- c("Form") 
      val <- data.frame()
      for(i in 1:length(inputs[,1])){
      res <- car::deltaMethod(object=coefs, g=as.character(inputs$Form[i]), 
                      vcov.=vcMat, level = 0.95 )
      val <- rbind(val, res)
     }
     row.names(val) <- 1:length(inputs[,1])
     val <- cbind(inputs, val)
     # lisRes <- apply(func, 1, 
     #               function(x) car::deltaMethod(object=coefs, g=x, 
     #                    vcov.=vcMat, level = 0.95))
     # val <- plyr::ldply(lisRes)
     # val <- cbind(func, val)
   }else{
     # Una sola costante è possibile, anche se può assumere più valori
     #inputs <- expand.grid(func$form, const[,1])
     inputs <- apply(const, 2, function(col) expand.grid(func$form, col))
     out <- inputs[[1]]
     if(length(const[1,]) > 1){
     for(i in 2:length(const[1,])){
       out <- cbind(out, inputs[[i]][,2]) 
     } }
     inputs <- out
     colnames(inputs) <- c("Form", colnames(const)) 
     val <- data.frame()
     for(i in 1:length(inputs[,1])){
     valCost <- inputs[i, 2:length(inputs[1,])]
     names(valCost) <- colnames(inputs)[2:length(inputs[1,])]
     # print(valCost)
     res <- car::deltaMethod(object=coefs, g=as.character(inputs$Form[i]), 
                      vcov.=vcMat, level = 0.95, 
                      constants = valCost)
     val <- rbind(val, res)
     }
     row.names(val) <- 1:length(inputs[,1])
     val <- cbind(inputs, val)
     # val
     # Corrections made on 4/11/2020 - Andrea
     # lisRes <- apply(const, 1, function(y) apply(func, 1, 
     #               function(x) car::deltaMethod(object=coefs, g=x, 
     #                    vcov.=vcMat, level = 0.95, constants = y )))
     # 
     # val <- plyr::ldply(do.call(rbind, lisRes))
     # funList <- func %>% slice(rep(1:n(), each = length(const[,1])))
     # func
     # constList <- const %>% slice(rep(1:n(), length(func[,1])))
     # val <- cbind(funList, constList, val)
   }
  val
  lenVal <- length(val[1,]) 
  retDF <- val[, c(-(lenVal-1), -lenVal)]
  retDF$"t-value" <- abs(retDF$Estimate/retDF$SE)
  resDF <- ifelse(is.null(dfr), Inf, dfr)
  retDF$"p-value" <- 2 * pt(retDF$"t-value", resDF, lower.tail = F)
  if(resDF == Inf) colnames(retDF)[length(colnames(retDF)) - 1] <- "Z-value"
  return(retDF)
}




