CVA<-function(dataset,groups){

##dataset is a multidimensional matrix of observations
##groups is a vector coding for groupings
## This part performs a standardisation of data	
  	numdati<-length(dataset[,1])
  	numcolonne<-length(dataset[1,])
  	numclasses<-length(levels(groups))
  	if (numcolonne>numclasses) numcanonical=numclasses-1 else numcanonical=numcolonne	
  	dataset.medie<-matrix(c(as.vector(manova(dataset~1)$coefficients)),numdati,numcolonne,byrow=TRUE)
  	dataset.centrata<-dataset-dataset.medie
  	prodotto<-t(dataset.centrata)%*%dataset.centrata
  	i<-array(c(1:numcolonne,1:numcolonne),dim=c(numcolonne,2))
  	dev<-prodotto[i]
  	devst<-sqrt(dev/(numdati-1))
  	dataset.devst<-matrix(c(devst),numdati,numcolonne,byrow=TRUE)
  	dataset.standardizzata<-dataset.centrata/dataset.devst
   #dataset.standardizzata <- dataset
## This part calculates eigenvalues and eigenvectors
  	maovst<-manova(dataset.standardizzata~groups)
  	fitted<-maovst$fitted
  	residui<-maovst$residuals
  	T<-t(dataset.standardizzata)%*%dataset.standardizzata
  	B<-t(fitted)%*%fitted
  	W<-t(residui)%*%residui
  	WB<-solve(W)%*%B
  	A<-as.numeric(eigen(WB)$values[1:numcanonical])
  	V1<-as.numeric(eigen(WB)$vectors[,1:numcanonical])
  	V1<-matrix(V1,numcolonne,numcanonical)
  	canvarnames<-c(1:numcanonical)
	canvarnames<-paste("CAN", canvarnames,sep="")
  
## eigenvalues are scaled to obtain standardised canonical coefficients
	VARCAN1<-dataset.standardizzata%*%V1
	maovst2<-manova(VARCAN1~groups)
	varcovar<-t(maovst2$residuals)%*%maovst2$residuals
	i<-array(c(1:numcanonical,1:numcanonical),dim=c(numcanonical,2))
			dev1<-varcovar[i]
			devst1<-sqrt(dev1/(numdati-length(levels(groups))))
			scaling<-1/devst1
			identita<-matrix(c(0),numcanonical,numcanonical)
			identita[i]<-scaling
		coefst<-V1%*%identita
	dimnames(coefst)<-list(dimnames(dataset)[[2]],canvarnames)

  ## calculation of raw canonical coefficients
	dev2<-t(dataset.centrata)%*%dataset.centrata
	i<-array(c(1:numcolonne,1:numcolonne),dim=c(numcolonne,2))
	dev2<-dev2[i]
	devst2<-sqrt(dev2/(numdati))
	devst2<-1/devst2
	identita<-matrix(c(0),numcolonne,numcolonne)
	identita[i]<-devst2
	raw<-identita%*%coefst
	dimnames(raw)<-list(dimnames(dataset)[[2]],canvarnames)

 ## canonical analysis
	prop<-A/sum(A)
	cancor<-sqrt(A/(1+A))
	cancor2<-cancor^2

 ## scores of centroids
	medie<-aov(dataset.standardizzata~groups-1)$coefficients
  	scores<-medie%*%coefst

  ## canonical structure
	VARCAN<-dataset.standardizzata%*%coefst
	total<-cor(dataset,VARCAN)

	medie<-aov(dataset.centrata~groups-1)$fitted
	VARCANmedie<-medie%*%raw
	between<-cor(medie,VARCANmedie)

	medie<-aov(dataset~groups-1)$residuals
	VARCANresidui<-medie%*%raw
	within<-cor(medie,VARCANresidui)

  ## RESULTS
	list("eigenvalues"=A, 
	     "eigenvectors"=V1, 
	     "Proportion" = prop, 
	     "Canonical correlation"=cancor, 
	     "Squared canonical correlation"=cancor2, 
	     "Standardised canonical coefficients"=coefst, 
	     "Raw canonical coefficients"=raw,
	     "Scores of centroids"=scores, 
	     "Total canonical structure"=total, 
	     "Between canonical structure"=between, 
	     "Within canonical structure" = within)
	
}