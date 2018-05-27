standardise<-function(dataset){

  ##dataset is a multidimensional matrix of observations	
  	numdati<-length(dataset[,1])
  	numcolonne<-length(dataset[1,])
  	dataset.medie<-matrix(c(as.vector(manova(dataset~1)$coefficients)),numdati,numcolonne,byrow=TRUE)
  	dataset.centrata<-dataset-dataset.medie
  	prodotto<-t(dataset.centrata)%*%dataset.centrata
  	i<-array(c(1:numcolonne,1:numcolonne),dim=c(numcolonne,2))
  	dev<-prodotto[i]
  	devst<-sqrt(dev/(numdati-1))
  	dataset.devst<-matrix(c(devst),numdati,numcolonne,byrow=TRUE)
  	dataset.standardizzata<-dataset.centrata/dataset.devst
  	dataset.standardizzata	
}