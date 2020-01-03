gnlht <- function(obj, func,  const = NULL){
   # obj <- modNlin
   coefs <- coef(obj)
   vcovC <- vcov(obj)
   # const <- constList; func <- funList
   # print(is.null(const))
   if(is.null(const) == T){
     lisRes <- apply(func, 1, 
                   function(x) car::deltaMethod(object=coefs, g=x, 
                        vcov.=vcovC, level = 0.95))
     val <- ldply(lisRes)
     val <- cbind(func, val)
   }else{
     lisRes <- apply(const, 1, function(y) apply(func, 1, 
                   function(x) car::deltaMethod(object=coefs, g=x, 
                        vcov.=vcovC, level = 0.95, constants = y )))
     val <- ldply(do.call(rbind, lisRes))
     funList <- func %>% slice(rep(1:n(), each = length(const[,1])))
     constList <- const %>% slice(rep(1:n(), length(func[,1])))
     val <- cbind(funList, constList, val)
   }
  return(val)
}

