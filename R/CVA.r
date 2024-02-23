CVA <- function(dataset, groups, scale = TRUE, constraint = 3){
## Performs canonical variate analysis (date: 15/11/19)
## dataset is a multidimensional data.frame of observations
## groups is a vector coding for groupings
## This part performs a standardisation of data	
    groups <- as.factor(groups)
    dataset <- as.matrix(dataset)
  	numdati <- length(dataset[,1])
  	numcolonne <- length(dataset[1,])
  	numclasses <- length(levels(groups))
  	if (numcolonne > numclasses) numcanonical = numclasses - 1 else numcanonical=numcolonne	
  	canvarnames<-c(1:numcanonical)
	  canvarnames<-paste("CAN", canvarnames,sep="")
  
  	# Centratura e calcolo coefficienti canonici raw
  	# dataset.centrata <- scale(dataset, center = TRUE, scale = F)
  	# Tc <- t(dataset.centrata) %*% dataset.centrata
  	# i <- array(c(1:numcolonne,1:numcolonne),dim=c(numcolonne,2))
  	# dev <- Tc[i]
  	# devst <- sqrt(dev/(numdati-1))
  	# dataset.devst <- matrix(c(devst),numdati,numcolonne,byrow=TRUE)
  	# 
  	# maovst <- manova(dataset.centrata ~ groups)
  	# fitted <- maovst$fitted
  	# residui <- maovst$residuals
  	# Bc <- t(fitted)%*%fitted
  	# Wc <- t(residui)%*%residui
  	# WBc <- solve(Wc) %*% Bc
  	# V1c <- as.numeric(eigen(WBc)$vectors[,1:numcanonical])
  	# V1c <- matrix(V1c, numcolonne, numcanonical)
  	# scal <- sqrt( diag( t(V1c) %*% (Wc/(numdati - numclasses)) %*% V1c ) )
  	# raw <- t( apply(V1c, 1, function(x) x / scal) )
  	# if(nrow(raw) < 2) raw <- t(raw)
  	# dimnames(raw) <- list(dimnames(dataset)[[2]], canvarnames)
  	# # raw
  	
  	# Eventuale tandardizzazione e calcolo dei coefficienti canonici
  	dataset.input <- scale(dataset, center = TRUE, scale = scale)
  	# print(dataset.input)
  	
	  maovst <- manova(dataset.input ~ groups)
  	fitted <- maovst$fitted
  	residui<-maovst$residuals
  	Ts <- t(dataset.input)%*%dataset.input
  	B <- t(fitted)%*%fitted
  	W <- t(residui)%*%residui
  	WB <- solve(W) %*% B
  	
  	## This part calculates eigenvalues and eigenvectors
  	A <- as.numeric(eigen(WB)$values[1:numcanonical])
  	V1 <- as.numeric(eigen(WB)$vectors[,1:numcanonical])
  	V1 <- matrix(V1, numcolonne, numcanonical)
  	
	  VARCAN1 <- dataset.input %*% V1
	  devstVarCan <- apply(VARCAN1, 2, sd)
	  if(length(VARCAN1[1,]) > 1) {
	       maovst2 <- manova(VARCAN1 ~ groups)
	  } else {
	       maovst2 <- lm(VARCAN1 ~ groups)
	  }
	  varcovar <- t(maovst2$residuals)%*%maovst2$residuals
	  i <- array(c(1:numcanonical,1:numcanonical),dim=c(numcanonical,2))
		dev1 <- varcovar[i]
		devst1 <- sqrt(dev1/(numdati - length(levels(groups))))
		if(constraint == 1){
		  # unit total sd
		  scaling <- 1/devstVarCan
		} else if(constraint == 2) {
		  # unit total norm
		  scaling <- 1/(devstVarCan * sqrt(numdati - 1))
		} else if(constraint == 3) {
		  # unit within sd
		  scaling <- 1/devst1
		} else {
		  # unit within norm
		  scaling <- 1/(devst1 * sqrt(numdati - length(levels(groups))))
		}
		# print(scaling)
		identita <- matrix(c(0),numcanonical,numcanonical)
		identita[i] <- scaling
		# print(V1)
		coefst <- V1 %*% identita
	  dimnames(coefst) <- list(dimnames(dataset)[[2]], canvarnames)
  	# coefst

## calculation of raw canonical coefficients
	  # dev2 <- t(dataset.centrata)%*%dataset.centrata
	  # i <- array(c(1:numcolonne,1:numcolonne),dim=c(numcolonne,2))
	  # dev2 <- dev2[i]
	  # devst2 <- sqrt(dev2/(numdati))
	  # devst2 <- 1/devst2
	  # identita <- matrix(c(0),numcolonne,numcolonne)
	  # identita[i]<-devst2
	  # raw <- identita%*%coefst
	  # dimnames(raw) <- list(dimnames(dataset)[[2]],canvarnames)

 ## canonical analysis
	   prop <- A/sum(A)
	   cancor <- sqrt(A/(1+A))
	   cancor2 <- cancor^2

 ## scores of centroids
	   medie <- aov(dataset.input ~ groups-1)$coefficients
     scores <- medie %*% coefst
     rownames(scores) <- gsub("groups", "", rownames(scores))
     
 ## Classification functions
     medie <- aov(as.matrix(dataset.input) ~ groups - 1)$coefficients
     Vcov <- W/(numdati - numclasses)
     constants <- diag(- 0.5 * medie %*% solve(Vcov) %*% t(medie) )
     class.fun <- cbind(constants, t(solve(Vcov) %*% t(medie)) )
      
     matDati <- as.matrix( cbind(rep(1,numdati), dataset.input))
     matCoef <- as.matrix( class.fun )
     classVal <- matDati %*% t(matCoef)
     class <- levels(groups)[apply(classVal, 1, function(x) which.max(x))]
    
  ## canonical scores
     VARCAN <- dataset.input%*%coefst
     row.names(VARCAN) <- groups
	   
	   
  ## canonical structure
	   total <- cor(dataset, VARCAN)

	   medie <- aov(dataset.input~groups-1)$fitted
	   VARCANmedie <- medie %*% coefst
	   between <- cor(medie,VARCANmedie)

	   medie <- aov(as.matrix(dataset.input) ~ groups-1)$residuals
	   VARCANresidui <- medie %*% coefst
	   within <- cor(medie, VARCANresidui)

  ## RESULTS
	list(
	     # MANOVA
	     "TOT" = Ts,
	     "B" = B,
	     "W" = W,
	     "B/W" = WB,
	     # SPECTRAL DECOMPOSITION
	     "eigenvalues" = A, 
	     "eigenvectors" = V1,
	     "proportion" = prop,
	     # CANONICAL ANALYSES
	     "correlation" = cancor, 
	     "squared.canonical.correlation" = cancor2, 
	     "coefficients"=coefst, 
	     # "Raw.coefficients"=raw,
	     "scores" = VARCAN,
	     "centroids"=scores,
	     # Canonical structure
	     "total.structure"=total, 
	     "between.structure"=between, 
	     "within.structure" = within,
	     "Class.fun" = class.fun,
	     "Class.val" = classVal,
	     "Class" = class
	     )
	
}
