gnlht <- function(obj, func,  const = NULL, vcov. = vcov, parameterNames=names(coef(obj))){
   # obj <- modNlin
   # func <- c(~k1 - k2, ~k1 - k3, ~k1 - k4)
   # Da qui
   temp <- lapply(func, function(x) as.character(as.expression(x[[length(x)]])))
   func <- data.frame(form=unlist(temp))
   coefs <- coef(obj)
   names(coefs) <- parameterNames
   #print(str(vcov.))
   if (is.function(vcov.)){
            vcMat <- vcov.(obj)
        } else {vcMat <- vcov.}
   
   # const <- constList; func <- funList
   # print(is.null(const))
   if(is.null(const) == T){
     lisRes <- apply(func, 1, 
                   function(x) car::deltaMethod(object=coefs, g=x, 
                        vcov.=vcMat, level = 0.95))
     val <- plyr::ldply(lisRes)
     val <- cbind(func, val)
   }else{
     lisRes <- apply(const, 1, function(y) apply(func, 1, 
                   function(x) car::deltaMethod(object=coefs, g=x, 
                        vcov.=vcMat, level = 0.95, constants = y )))
     val <- ldply(do.call(rbind, lisRes))
     funList <- func %>% slice(rep(1:n(), each = length(const[,1])))
     constList <- const %>% slice(rep(1:n(), length(func[,1])))
     val <- cbind(funList, constList, val)
   }
  retDF <- val[,1:3]
  retDF$"t-value" <- retDF$Estimate/retDF$SE
  if(class(obj) == "nls") resDF <- summary(obj)$df[2]
  else if(class(obj) == "drc") resDF <- summary(obj)$df
  else resDF <- Inf
  print(resDF)
  retDF$"p-value" <- 2 * pt(retDF$"t-value", resDF, lower.tail = F)
  return(retDF)
}

