GGE_old<-function(variety,envir,block,yield,PC=2){

  ## GENOTYPE PLUS GENOTYPE X INTERACTION, AS APPLIED TO AGRICULTURE VARIETY TRIALS
  ## REFERENCE: Yan et al., 2000. Cultivar evaluation and mega-environment investigation
  ## based on the GGE biplot
  ## Update: 30/10/2007

  ## INPUTS:
  ## variety is a vector of variety codes (factor)
  ## envir is a vector of environment codes (factor)
  ## block is a vector of block codes (factor)
  ## yield is a vector of yields to be analysed (numerical)
  ## PC is the number of PC to be considered (default is 2)

  ## 1 - Descriptive statistics
  overall.mean<-mean(yield)
  envir.mean<-tapply(yield,envir,mean)
  var.mean<-tapply(yield,variety,mean)  
  int.mean<-tapply(yield,list(variety,envir),mean)
  envir.num<-length(envir.mean)
  var.num<-length(var.mean)
  
  ## 2 - ANOVA (additive model) and table of GE effects
  variety<-factor(variety)
  envir<-factor(envir)
  block<-factor(block)  
  add.anova<-aov(yield~envir+block%in%envir+variety+envir:variety)
  add.anova.residual.SS<-deviance(add.anova)
  add.anova.residual.DF<-add.anova$df.residual
  add.anova.residual.MS<-add.anova.residual.SS/add.anova.residual.DF
  anova.table<-summary(add.anova)
  row.names(anova.table[[1]])<-c("Environments","Genotypes","Blocks(Environments)","Environments x Genotypes", "Residuals")
  ge.anova<-aov(yield~envir+envir:variety)
  ge.eff<-model.tables(ge.anova,type="effects",cterms="envir:variety")$tables$"envir:variety"
  
  ## 3 - SVD
  dec <- svd(ge.eff, nu = PC, nv = PC)
  if (PC > 1){ 
    D <- diag(dec$d[1:PC]) 
  } else {
    D <- dec$d[1:PC]
  }
  E<-dec$u%*%sqrt(D)
  G<-dec$v%*%sqrt(D)
  Ecolnumb<-c(1:PC)
  Ecolnames<-paste("PC",Ecolnumb,sep="")
  dimnames(E)<-list(levels(envir),Ecolnames)
  dimnames(G)<-list(levels(variety),Ecolnames)
 
  ## 4 - Singular values and significance of PCs  
  svalues<-dec$d
  PC.n<-c(1:length(svalues))
  PC.SS<-svalues^2
  percSS<-PC.SS/sum(PC.SS)*100
  GGE.table<-data.frame("PC"=PC.n,"Singular_value"=svalues,"PC_SS"=PC.SS, "Perc_of_Total_SS"=percSS)
  numblock<-length(levels(block))
  GGE.SS<-(t(as.vector(ge.eff))%*%as.vector(ge.eff))*numblock
  
  ## 5 - GGE Biplot
  plot(1,type='n',xlim=range(c(E[,1],G[,1])),ylim=range(c(E[,2],G[,2])),xlab="PC 1",ylab="PC 2")
  points(E[,1],E[,2],"n",col="red",lwd=5)
  text(E[,1],E[,2],labels=row.names(E),adj=c(0.5,0.5),col="red")
  points(G[,1],G[,2], "n", col="blue",lwd=5)
  text(G[,1],G[,2],labels=row.names(G),adj=c(0.5,0.5),col="blue")

  ## 6 - Results
  result <- list("Genotype_means"=var.mean, "Environment_means"=envir.mean, "Interaction_means"=int.mean, 
  "GE_effect"=ge.eff, "Additive_ANOVA"=anova.table, "GGE_Sum of Squares"=GGE.SS, GGE_summary = GGE.table, Environment_scores=E, "Genotype_scores"=G)
  class(result) <- "GGEobject"
  return(result)
}

GGEmeans <- function(variety, envir, yield, PC=2, scale=0.5){

  ## GENOTYPE PLUS GENOTYPE X INTERACTION, AS APPLIED TO AGRICULTURE VARIETY TRIALS
  ## REFERENCE: Yan et al., 2000. Cultivar evaluation and mega-environment investigation
  ## based on the GGE biplot
  ## Update: 30/11/2008

  ## INPUTS:
  ## variety is a vector of variety codes (factor)
  ## envir is a vector of environment codes (factor)
  ## yield is a vector of yields to be analysed (numerical)
  ## PC is the number of PC to be considered (default is 2)
  ## scale is the partitioning of singular values

  ge.anova <- aov(yield ~ envir + envir:variety)
  ge.anova2 <- aov(yield ~ envir + variety + envir:variety)
  ge.eff<-model.tables(ge.anova,type="effects",cterms="envir:variety")$tables$"envir:variety"
  means.yield<-model.tables(ge.anova, type="means", cterms="envir:variety")$tables$"envir:variety"
  means.env <- model.tables(ge.anova2, type="means", cterms="envir")$tables$"envir"
  means.gen <- model.tables(ge.anova2, type="means", cterms="variety")$tables$"variety"  
  
  ## 3 - SVD
  ge.eff <- t(ge.eff)
  dec <- svd(ge.eff, nu = PC, nv = PC)
    if (PC > 1){ 
    D <- diag(dec$d[1:PC]^scale)
    Dbis <- diag(dec$d[1:PC]^(1-scale)) 
  } else {
    D <- dec$d[1:PC]^scale
    Dbis <- dec$d[1:PC]^(1-scale)
    
  }
  
  G<-dec$u%*%D
  E<-dec$v%*%Dbis
  
  Ecolnumb<-c(1:PC)
  Ecolnames<-paste("PC",Ecolnumb,sep="")
   dimnames(E)<-list(levels(envir),Ecolnames)
   dimnames(G)<-list(levels(variety),Ecolnames)
 
  ## 4 - Singular values and significance of PCs  
  svalues<-dec$d
  PC.n<-c(1:length(svalues))
  PC.SS<-svalues^2
  percSS<-PC.SS/sum(PC.SS)*100
  GGE.table<-data.frame("PC"=PC.n,"Singular_value"=svalues,"PC_SS"=PC.SS, "Perc_of_Total_SS"=percSS)
  #numblock<-length(levels(block))
  #GGE.SS<-(t(as.vector(ge.eff))%*%as.vector(ge.eff))*numblock
  GGE.SS<-(t(as.vector(ge.eff))%*%as.vector(ge.eff))
  
  ## 6 - Results
  result <- list(
  Yield_means=t(means.yield), GE_effect=ge.eff, "GGE_Sum of Squares"=GGE.SS, 
  GGE_summary = GGE.table, Environment_scores=E, "Genotype_scores"=G, "Genotype_means"=means.gen, "Environment_means"=means.env)
  class(result) <- "GGEobject"
  cat(paste("Result of GGE Analysis", "\n", "\n"))
	cat(paste("Environment Scores", "\n"))
	print(E)
	cat(paste("\n","Genotype Scores", "\n"))
	print(G)
	return(invisible(result))
}

GGE<-function(variety, envir, block, yield, PC=2, scale=0.5){

  ## GENOTYPE PLUS GENOTYPE X INTERACTION, AS APPLIED TO AGRICULTURE VARIETY TRIALS
  ## REFERENCE: Yan et al., 2000. Cultivar evaluation and mega-environment investigation
  ## based on the GGE biplot
  ## Update: 30/10/2007

  ## INPUTS:
  ## variety is a vector of variety codes (factor)
  ## envir is a vector of environment codes (factor)
  ## block is a vector of block codes (factor)
  ## yield is a vector of yields to be analysed (numerical)
  ## PC is the number of PC to be considered (default is 2)
  ## scale is the partitioning of singular values

  ## 1 - Descriptive statistics
  overall.mean<-mean(yield)
  envir.mean<-tapply(yield,envir,mean)
  var.mean<-tapply(yield,variety,mean)  
  int.mean<-tapply(yield,list(variety,envir),mean)
  envir.num<-length(envir.mean)
  var.num<-length(var.mean)
  
  ## 2 - ANOVA (additive model) and table of GE effects
  variety<-factor(variety)
  envir<-factor(envir)
  block<-factor(block)  
  add.anova<-aov(yield~envir+block%in%envir+variety+envir:variety)
  add.anova.residual.SS<-deviance(add.anova)
  add.anova.residual.DF<-add.anova$df.residual
  add.anova.residual.MS<-add.anova.residual.SS/add.anova.residual.DF
  anova.table<-summary(add.anova)
  row.names(anova.table[[1]])<-c("Environments","Genotypes","Blocks(Environments)","Environments x Genotypes", "Residuals")
  ge.anova<-aov(yield~envir+envir:variety)
  ge.eff<-model.tables(ge.anova,type="effects",cterms="envir:variety")$tables$"envir:variety"
  ge.eff <- t(ge.eff)
  
  ## 3 - SVD
  dec <- svd(ge.eff, nu = PC, nv = PC)
  if (PC > 1){ 
    D <- diag(dec$d[1:PC]^scale)
    Dbis <- diag(dec$d[1:PC]^(1-scale)) 
  } else {
    D <- dec$d[1:PC]^scale
    Dbis <- dec$d[1:PC]^(1-scale)
    
  }
  
  G<-dec$u%*%D
  E<-dec$v%*%Dbis
  
  Ecolnumb<-c(1:PC)
  Ecolnames<-paste("PC",Ecolnumb,sep="")
  dimnames(E)<-list(levels(envir),Ecolnames)
  dimnames(G)<-list(levels(variety),Ecolnames)
 
  ## 4 - Singular values and significance of PCs  
  svalues<-dec$d
  PC.n<-c(1:length(svalues))
  PC.SS<-svalues^2
  percSS<-PC.SS/sum(PC.SS)*100
  GGE.table<-data.frame("PC"=PC.n,"Singular_value"=svalues,"PC_SS"=PC.SS, "Perc_of_Total_SS"=percSS)
  numblock<-length(levels(block))
  GGE.SS<-(t(as.vector(ge.eff))%*%as.vector(ge.eff))*numblock
  
  ## 5 - Results
  result <- list(Genotype_means=var.mean, Environment_means=envir.mean, Interaction_means=int.mean, 
  GE_effect=ge.eff, Additive_ANOVA=anova.table, "GGE_Sum of Squares"=GGE.SS, 
  GGE_summary = GGE.table, Environment_scores=E, "Genotype_scores"=G)
  class(result) <- "GGEobject"
    cat(paste("Result of GGE Analysis", "\n", "\n"))
	cat(paste("Environment Scores", "\n"))
	print(E)
	cat(paste("\n","Genotype Scores", "\n"))
	print(G)
	return(invisible(result))
}
