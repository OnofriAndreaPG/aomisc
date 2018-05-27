TukeyNA<-function(dati, tesi, blocco){
	## dati, tesi and blocco are vectors of respectively data, treatment codes and block codes)
	tesi<-factor(tesi)
	blocco<-factor(blocco)
	model<-lm(dati~tesi+blocco)
	model.analisi<-(anova(model))

	model.tesi<-lm(dati~tesi-1)
	model.blocco<-lm(dati~blocco-1)
	d1<-model.blocco$fitted-mean(dati)
	d2<-model.tesi$fitted-mean(dati)
	N<-sum(d1*d2*dati)
		D<-(sum(d1^2)/length(levels(tesi)))*(sum(d2^2)/length(levels(blocco)))

		NonAdditivity<-N^2/D
		F<-NonAdditivity/model.analisi$Mean[3]
		DF<-model.analisi$Df[3]-1
		remainder<-model.analisi$Sum[3]-NonAdditivity
	
		remainderMS<-remainder/DF

		F<-NonAdditivity/remainderMS

		probF<-1-pf(F,1,anova(model)$Df[3]-1)
		ris.1<-c(1,7)
		ris.2<-c(NonAdditivity,remainder)
		ris.3<-c(NonAdditivity,remainderMS)
		ris.4<-c(round(F,2),"")
		ris.5<-c(round(probF,7),"")
		B<-N/D
		p<-1-B*mean(dati)
		risultati<-data.frame(df=ris.1,"Sum of Squares"=ris.2,"Mean Squares"=ris.3,F=ris.4,Prob.=ris.5)
		row.names(risultati)<-c("Non additivity", "Residual")
		list("Test of Additivity"=risultati,"Suggested value for power transformation of data" =p)

}