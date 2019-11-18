CVA<-function(dataset, groups){
## Performs linear discriminant analysis (date: 15/11/19)
## dataset is a multidimensional data.frame of observations
## groups is a vector coding for groupings
## This part performs a standardisation of data	
  	numdati <- length(dataset[,1])
  	numcolonne <- length(dataset[1,])
  	numclasses <- length(levels(groups))
  	if (numcolonne > numclasses) numcanonical = numclasses - 1 else numcanonical=numcolonne	
  	canvarnames<-c(1:numcanonical)
	  canvarnames<-paste("CAN", canvarnames,sep="")
  
  	#Centratura e calcolo coefficienti canonici raw
  	dataset.centrata <- scale(dataset, scale = F)
  	Tc <- t(dataset.centrata)%*%dataset.centrata
  	i <- array(c(1:numcolonne,1:numcolonne),dim=c(numcolonne,2))
  	dev <- Tc[i]
  	devst <- sqrt(dev/(numdati-1))
  	dataset.devst <- matrix(c(devst),numdati,numcolonne,byrow=TRUE)
  	
  	maovst <- manova(dataset.centrata ~ groups)
  	fitted <- maovst$fitted
  	residui <- maovst$residuals
  	Bc <- t(fitted)%*%fitted
  	Wc <- t(residui)%*%residui
  	WBc <- solve(Wc) %*% Bc
  	V1c <- as.numeric(eigen(WBc)$vectors[,1:numcanonical])
  	V1c <- matrix(V1c, numcolonne, numcanonical)
  	scal <- sqrt( diag( t(V1c) %*% (Wc/(numdati - numclasses)) %*% V1c ) )
  	raw <- t( apply(V1c, 1, function(x) x / scal) )
  	if(nrow(raw) < 2) raw <- t(raw)
  	dimnames(raw) <- list(dimnames(dataset)[[2]], canvarnames)
  	raw
  	
  	#Standardizzazione e calcolo coefficienti standardizzati
  	dataset.standardizzata <- scale(dataset, scale = T)
    maovst <- manova(dataset.standardizzata ~ groups)
  	fitted <- maovst$fitted
  	residui<-maovst$residuals
  	Ts <- t(dataset.standardizzata)%*%dataset.standardizzata
  	B <- t(fitted)%*%fitted
  	W <- t(residui)%*%residui
  	WB <- solve(W)%*%B
  	
  	## This part calculates eigenvalues and eigenvectors
  	A <- as.numeric(eigen(WB)$values[1:numcanonical])
  	V1 <- as.numeric(eigen(WB)$vectors[,1:numcanonical])
  	V1 <- matrix(V1, numcolonne, numcanonical)
	  VARCAN1 <- dataset.standardizzata %*% V1
	  if(length(VARCAN1[1,]) > 1) {
	       maovst2 <- manova(VARCAN1 ~ groups)
	  } else {
	       maovst2 <- lm(VARCAN1 ~ groups)
	  }
	  varcovar <- t(maovst2$residuals)%*%maovst2$residuals
	  i <- array(c(1:numcanonical,1:numcanonical),dim=c(numcanonical,2))
		dev1 <- varcovar[i]
		devst1 <- sqrt(dev1/(numdati-length(levels(groups))))
		scaling <- 1/devst1
		identita <- matrix(c(0),numcanonical,numcanonical)
		identita[i] <- scaling
		coefst <- V1 %*% identita
	  dimnames(coefst) <- list(dimnames(dataset)[[2]], canvarnames)
  	coefst

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
	   cancor<-sqrt(A/(1+A))
	   cancor2<-cancor^2

 ## scores of centroids (creati con gli standardised coefficients)
	   medie <- aov(dataset.standardizzata~groups-1)$coefficients
     scores <- medie%*%coefst

 ## Classification functions
     medie <- aov(as.matrix(dataset) ~ groups - 1)$coefficients
     Vcov <- Wc/(numdati - numclasses)
     constants <- diag(- 0.5 * medie %*% solve(Vcov) %*% t(medie) )
     class.fun <- cbind(constants, t(solve(Vcov) %*% t(medie)) )
  
     matDati <- as.matrix( cbind(rep(1,numdati), dataset))
     matCoef <- as.matrix( class.fun )
     classVal <- matDati %*% t(matCoef)
     class <- levels(groups)[apply(classVal, 1, function(x) which.max(x))]
    

  ## canonical structure
	   VARCAN <- dataset.standardizzata%*%coefst
	   total <- cor(dataset, VARCAN)

	   medie <- aov(dataset.centrata~groups-1)$fitted
	   VARCANmedie <- medie %*% raw
	   between <- cor(medie,VARCANmedie)

	   medie <- aov(as.matrix(dataset) ~ groups-1)$residuals
	   VARCANresidui <- medie%*%raw
	   within <- cor(medie,VARCANresidui)

  ## RESULTS
	list("eigenvalues" = A, 
	     "eigenvectors.c" = V1c,
	     "eigenvectors.s" = V1,
	     "Proportion" = prop, 
	     "Canonical.correlation" = cancor, 
	     "Squared.canonical.correlation" = cancor2, 
	     "Standardised.coefficients"=coefst, 
	     "Raw.coefficients"=raw,
	     "Canonical.scores" = VARCAN,
	     "Centroids.scores"=scores, 
	     "Total.structure"=total, 
	     "Between.structure"=between, 
	     "Within.structure" = within,
	     "Class.fun" = class.fun,
	     "Class.val" = classVal,
	     "Class" = class)
	
}
