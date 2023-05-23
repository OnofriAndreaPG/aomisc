compCoefs <- function (object, parameterNames, operator = "-", 
                       display = "pairwise", 
                       adjust = "holm", reversed = FALSE, level = 0.05,
                       Letters = c(letters, LETTERS)) {
  # Make pairwise comparisons based on a 'vector of means'nls' object.
  # It is shaped after the drc::compParm function and depends only on 
  # multcompView. Adjustments for multiplicity are possible, using the
  # p.adjust method. Holm is the default
  # Updated on 19/05/2023
  
  if(!inherits(object, "nls")) stop("This function only works on 'nls' objects")
  # object <- modNlin
  # parameterNames <- pn
  
  # Retrieve model information
  nams <- names(object$m$getPars())
  selectedCoefs <- nams %in% parameterNames
  lenPV <- length(selectedCoefs[selectedCoefs==T])
  if (lenPV < 2) {
    stop("No parameters to compare were found")
  }
  
  strVal <- factor(parameterNames, levels = as.character(parameterNames))
  parm <- coef(object)[selectedCoefs]
  varMat <- vcov(object)[selectedCoefs, selectedCoefs]
  df.residual <- df.residual(object)
  
  # Sort the vectors/matrices (so that letters are in 'traditional' order)
  tmp <- data.frame(parm=parm, strVal=strVal)
  if(!reversed) varMat <- varMat[order(parm), order(parm)] else varMat <- varMat[order(parm, decreasing = TRUE), order(parm, decreasing = TRUE) ]
  if(!reversed) tmp <- tmp[order(parm), ] else tmp <- tmp[order(parm, decreasing = TRUE), ]
  parm <- tmp$parm; presentVec <- tmp$strVal # for better re-using the code in compParm()
  
  if (identical(operator, "/")) {
    hypVal <- 1
    fct <- function(ind) {
      parm[ind[1]]/parm[ind[2]]
    }
    dfct <- function(ind) {
      transVec <- c(1/parm[ind[2]], -parm[ind[1]]/(parm[ind[2]]^2))
      sqrt(transVec %*% varMat[ind, ind] %*% transVec)
    }
  }
  if (identical(operator, "-")) {
    hypVal <- 0
    fct <- function(ind) {
      parm[ind[1]] - parm[ind[2]]
    }
    transVec <- c(1, -1)
    dfct <- function(ind) {
      sqrt(transVec %*% varMat[ind, ind] %*% transVec)
    }
  }
  lenRV <- lenPV * (lenPV - 1)/2
  cpMat <- matrix(0, lenRV, 4)
  compParm <- rep("", lenRV)
  degfree <- ifelse(is.na(df.residual) == T, Inf, df.residual)
  pFct <- function(x) {
    pt(x, degfree)
  }
  strParm <- presentVec
  k <- 1
  for (i in 1:lenPV) {
    for (j in 1:lenPV) {
      if (j <= i) {
        next
      }
      ind <- which(presentVec %in% presentVec[c(i, j)])
      cpMat[k, 1] <- fct(ind)  # Diff
      cpMat[k, 2] <- dfct(ind) # SE
      tVal <- (cpMat[k, 1] - hypVal)/cpMat[k, 2]
      cpMat[k, 3] <- tVal
      cpMat[k, 4] <- pFct(-abs(tVal)) + (1 - pFct(abs(tVal)))
      compParm[k] <- paste(strParm[ind], 
                           collapse = operator)
      k <- k + 1
    }
  }
  cpMat <- data.frame(cpMat)
  row.names(cpMat) <- compParm
  colnames(cpMat) <- c("Estimate", "Std. Error", 
                       "t-value", "p-value")
  adjusted.P <- p.adjust(cpMat$`p-value`, method = as.character(adjust))
  cpMat$`p-value` <- as.vector(adjusted.P)
  
  if(display == "pairwise"){
    return(cpMat)
  } else {
    p.logic <- as.vector(adjusted.P)
    p.logic <- ifelse(p.logic > 0.05, FALSE, TRUE)
    names(p.logic) <- sub("/", "-", as.vector(row.names(cpMat)))
    tmp <- data.frame(parNam = presentVec, parm=parm)
    # if(reversed == T | reversed == F) {
   # Letters <- multcompView::multcompLetters3("parNam", "obj", p.logic,
   #                 data = tmp, reverse = F,
   #                 threshold = level)
   Letters <- multcompView::multcompLetters(p.logic, threshold = level,
                                       Letters = Letters, reversed = F)
   parMat <- data.frame(Value = parm, SE = sqrt(diag(varMat)))
   row.names(parMat) <- presentVec
   # parMat <- parMat[order(-parMat$Value),]
    # } 
    # else {
    #    # Letters <- multcompView::multcompLetters3("parNam", "obj", p.logic,
    #    #                 data = tmp, reverse = T,
    #    #                 threshold = level)
    #    Letters <- multcompView::multcompLetters(p.logic, threshold = level,
    #                                        Letters = Letters, reversed = T)
    #    parMat <- data.frame(Value = parm, SE = sqrt(diag(varMat)))
    #    row.names(parMat) <- parameterNames
    #    parMat <- parMat[order(-parMat$Value),]
    # }
   parMat$CLD = as.character(Letters$Letters)
  
  if(reversed == F) { parMat <- parMat[order(parMat$Value), ]
  } else {parMat <- parMat[order(-parMat$Value), ] }
   return(parMat)
  }
}
