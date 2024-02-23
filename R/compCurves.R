# Code to compare regression curves in a pairwise fashion
# The two functions are different: in the first one, the whole
# dataset is considered and a model is coded, where all the
# curves are different, except the two curves under comparison
# The second function works by considering subsets of the whole
# dataset. The first function appears to be more logic!
# Data: 22/01/2024

compCurves <- function(obj, adjusted = "none"){
  curveidVar <- obj$dataList$curveid
  if(inherits(obj, "drcte")) curveidVar <- obj$data[,length(obj$data) - 1]
  nlevs <- length(levels(curveidVar))
  comp <- combn(1:nlevs, 2)
  npairs <- length(comp[1,])
  pVals <- c(); RSS1 <- c(); RSS2 <- c(); DF <- c(); Fval <- c(); nam <- c()
  
  for(pair in 1:npairs) {
    # pair <- 1
    dummy1 <- as.factor(ifelse(curveidVar == levels(curveidVar)[comp[1,pair]] |
         curveidVar == levels(curveidVar)[comp[2,pair]],
         "D1", as.character(curveidVar)))
    obj2 <- update(obj, curveid = dummy1)
    if(inherits(obj, "drcte")) {
      aov <- anova(obj, obj2, details = FALSE, test = "Chisq")
      nam[pair] <- paste(levels(curveidVar)[comp[1,pair]], 
                 levels(curveidVar)[comp[2,pair]], sep = "-")
      RSS1[pair] <- aov$Loglik[1]
      RSS2[pair] <- aov$Loglik[2]
      DF[pair] <- aov$Df[2]
      Fval[pair] <- aov$`LR value`[2]
      pVals[pair] <- aov$`p value`[2]
    } else {
      aov <- anova(obj, obj2, details = FALSE)
      nam[pair] <- paste(levels(curveidVar)[comp[1,pair]], 
                 levels(curveidVar)[comp[2,pair]], sep = "-")
      RSS1[pair] <- aov$RSS[1]
      RSS2[pair] <- aov$RSS[2]
      DF[pair] <- aov$Df[2]
      Fval[pair] <- aov$`F value`[2]
      pVals[pair] <- aov$`p value`[2]
      }
  }
  pVals <- p.adjust(pVals, method = adjusted)
  dfRes <- data.frame(RSS1, RSS2, DF, "Test value"=Fval, "p-level" = pVals)
  row.names(dfRes) <- nam
  
  # Add letters
  p.logic <- ifelse(pVals > 0.05, FALSE, TRUE)
  names(p.logic) <- nam
  Letters <- multcompLetters(p.logic)
  dfLett <- data.frame("Curve" = levels(factor(curveidVar)), "Letters" = Letters$Letters) 
  return(list("Pairs" = dfRes, "Letters" = dfLett))
}

compCurves.2 <- function(obj, adjusted = "none"){
    df <- obj$data
    curveId <- levels(factor(df[,4]))
    el <- length(curveId)
    first <- c(); second <- c()
    for(i in 1:(el-1)) for(j in (i + 1):el){
      first <- c(first, i)
      second <- c(second, j)
    }
    cfr <- data.frame(first, second)
    pr <- apply(cfr, 1, function(x) {
      RedSub <- df[df[,4] == levels(df[,4])[x[1]] |
                   df[,4] == levels(df[,4])[x[2]], ]
      m1 <- update(obj, data = RedSub)
      m2 <- update(m1, curveid = NULL)
      tmp <- anova(m1, m2, details = F)
      pval <- tmp[2,5]
      RSS1 <- tmp[1,2]
      RSS2 <- tmp[2,2]
      DF <- tmp[2,3]
      Fval <- tmp[2,4]
      namC <- paste(x[1], x[2], sep = "-")
      data.frame(namC = namC, RSS1, RSS2, DF, Fval, pval = pval)
    })
    pr <- do.call(rbind, pr)
    pr[,6] <- p.adjust(pr[,6], method = adjusted)
    Probs <- pr[,6]
    
    # names(Probs) <- pr[,2]
    # Probs
    p.logic <- ifelse(Probs > 0.05, FALSE, TRUE)
    names(p.logic) <- pr[,1]
    Letters <- multcompLetters(p.logic)
    letRes <- data.frame("Curve" = levels(factor(df[,4])), "Letters" = Letters$Letters) 
    return(list("Pairs" = pr, "Letters" = letRes))
}