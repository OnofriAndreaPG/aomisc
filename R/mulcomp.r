mulcomp <- function (object, ...) UseMethod("mulcomp")

mulcomp.numeric <- function(fitted, nis, df, MSE, conf.level=.05, type="LSD", decreasing=TRUE){
  ## fitted is a vector of means
  ## nis is the number of replicates
  ## df is the residual df from the ANOVA table
  ## MSE is the mean squared error from the ANOVA table
  ## conf.level is the family-wise confidence level, defaults to .05
  ## type is the type of text
  ## decreasing is the type of sorting, default to TRUE
  num.means<-length(fitted)
  if(length(nis) == 1) nis <- rep(nis,num.means)  
  names(fitted)<-paste(c("M"),1:length(fitted),sep="")
  if(decreasing==TRUE){
  fitted<-sort(fitted, decreasing=TRUE)
  }
  else{
  fitted<-sort(fitted, decreasing=FALSE)
  }
  
  
  diffs<-matrix(0,num.means,num.means)
  for(i in 1:num.means){
    for(j in 1:num.means){
      diffs[i,j]<-fitted[j]-fitted[i]
    }
  }


  if (type=="LSD") {crit.diff<-LSD(nis, df, MSE, conf.level)}
  if (type=="NewmanKeuls"){crit.diff<-NewmanKeuls(nis, df, MSE, conf.level)} 
  if (type=="Duncan"){crit.diff<-Duncan(nis, df, MSE, conf.level)}
  if (type=="TukeyHSD"){crit.diff<-TukHSD(nis, df, MSE, conf.level)}
  if (type=="TukeyMRT"){crit.diff<-TukMRT(nis, df, MSE, conf.level)}
  if (type=="Scheffe"){crit.diff<-Scheffe(nis, df, MSE, conf.level)}
  if (type=="Bonferroni"){crit.diff<-Bonferroni(nis, df, MSE, conf.level)}
  if (type=="Sidak"){crit.diff<-Sidak(nis, df, MSE, conf.level)}
  
  sign<-matrix("",nrow=num.means,ncol=num.means)
  ext.lett<-c(letters,LETTERS)
  lettera<-1
  col<-1
  pos<-0
  posold<-0
  
 while(col<=num.means){
    if(abs(diffs[num.means,col])<=crit.diff[num.means,col]){
       for(i in col:num.means){
          sign[i,col]<-ext.lett[lettera]          
       }
      break                                                       
    }
    else{
      pos<-num.means
      while(abs(diffs[pos,col])>=crit.diff[pos,col]){        
        pos<-pos-1
        if (pos==col) break       
        }
      if (pos!=posold){
        for(i in col:pos){
          sign[i,col]<-ext.lett[lettera]
        }
        lettera<-lettera+1
        posold<-pos
      }
    
    }
   col<-col+1 
    
  }
  sign.lett<-sapply(1:num.means,function(nc)paste(sign[nc,],collapse=''))
  result<-data.frame(medie=fitted, test=sign.lett, row.names=names(fitted))
  cat(paste("Multiple Comparison testing", "\n"))
  print(result)
  cat(paste("\n"))
  report <- list(sem=sqrt(MSE/nis), sed=sqrt(2*MSE/nis), MCT=result, "Critical.Differences"=crit.diff[,1])
  invisible(report)  
}

mulcomp.aov <- function(x, conf.level=.05, type="LSD", decreasing=TRUE){
  ## x is an 'aov' object
  ## ns is the number of replicates
  ## df is the residual df from the ANOVA table
  ## MSE is the mean squared error from the ANOVA table
  ## conf.level is the family-wise confidence level, defaults to .05
  ## type is the type of text
  ## decreasing is the type of sorting, default to TRUE
  
  #Verifica se è aov object
  if (class(x)[1] != "aov"){  #& class(x)[1] != "aovlist"
  cat(paste("This function may be used only on 'aov' objects", "\n"))}
  
  #if so, then
  else{
  n.terms <- length(model.tables(x, type="means")$tables)
  report <- list(0)
  for(k in 2:n.terms){
  name.term <- names(model.tables(x, type="means")$n)[k-1]
  fitted  <- model.tables(x, type="means")$tables[[k]]
  num.means <- length(fitted)
  nis <- as.vector(model.tables(x, type="means")$n[[k-1]])
  if(length(nis) == 1) nis <- rep(nis, num.means)
  #Qui cambia a seconda che l'oggetto è aov o aovlist
  if (class(x)[1] == "aov" ){
    df <- x$df.residual
    MSE <- deviance(x)/df
    }
  else{
    df <- x[[k]]$df.residual
    MSE <- deviance(x[[k]])/df
    }
  
  #names(fitted)<-paste(c("M"),1:length(fitted),sep="")
  if(decreasing==TRUE){
  fitted<-sort(fitted, decreasing=TRUE)
  }
  else{
  fitted<-sort(fitted, decreasing=FALSE)
  }
   
  diffs<-matrix(0,num.means,num.means)
  for(i in 1:num.means){
    for(j in 1:num.means){
      diffs[i,j]<-fitted[j]-fitted[i]
    }
  }

  
  if (type=="LSD") {crit.diff<-LSD(nis, df, MSE, conf.level)}
  if (type=="NewmanKeuls"){crit.diff<-NewmanKeuls(nis, df, MSE, conf.level)} 
  if (type=="Duncan"){crit.diff<-Duncan(nis, df, MSE, conf.level)}
  if (type=="TukeyHSD"){crit.diff<-TukHSD(nis, df, MSE, conf.level)}
  if (type=="TukeyMRT"){crit.diff<-TukMRT(nis, df, MSE, conf.level)}
  if (type=="Scheffe"){crit.diff<-Scheffe(nis, df, MSE, conf.level)}
  if (type=="Bonferroni"){crit.diff<-Bonferroni(nis, df, MSE, conf.level)}
  if (type=="Sidak"){crit.diff<-Sidak(nis, df, MSE, conf.level)}
  
  sign<-matrix("",nrow=num.means,ncol=num.means)
  ext.lett<-c(letters,LETTERS)
  lettera<-1
  col<-1
  pos<-0
  posold<-0
  while(col<=num.means){
    
    if(abs(diffs[num.means,col])<=crit.diff[num.means,col]){
    
       for(i in col:num.means){
          sign[i,col]<-ext.lett[lettera]          
       }
      break                                                       
    }
    else{
      pos<-num.means
      while(abs(diffs[pos,col])>=crit.diff[pos,col]){        
        pos<-pos-1
        if (pos==col) break       
        }
      if (pos!=posold){
        for(i in col:pos){
          sign[i,col]<-ext.lett[lettera]
        }
        lettera<-lettera+1
        posold<-pos
      }
    
    }
   col<-col+1 
    
  }
  
  sign.lett<-sapply(1:num.means,function(nc)paste(sign[nc,],collapse=''))
  result<-data.frame(medie=fitted, test=sign.lett, row.names=names(fitted))  
  cat(paste("term: ", name.term, "\n"))
  print(result)
  cat(paste("\n"))
  #print("Critical.Differences"=crit.diff[,1])
  report[[k-1]] <- list(Effect=name.term, sem=sqrt(MSE/nis), sed=sqrt(2*MSE/nis), MCT=result, "Critical.Differences"=crit.diff[,1]) 
}
invisible(report)
}   
}

mulcomp.aovlist <- function(x, conf.level=.05, type="LSD", decreasing=TRUE){
  ## x is an 'aov' object
  ## ns is the number of replicates
  ## df is the residual df from the ANOVA table
  ## MSE is the mean squared error from the ANOVA table
  ## conf.level is the family-wise confidence level, defaults to .05
  ## type is the type of text
  ## decreasing is the type of sorting, default to TRUE
  
  #Verifica se è aov object
  if (class(x)[1] != "aovlist"){  #& class(x)[1] != "aovlist"
  cat(paste("This function may be used only on 'aov' objects", "\n"))}
  
  #if so, then
  else{
  n.terms <- length(model.tables(x, type="means")$tables)
  se.terms <- unlist(model.tables(x, type="effects", se=T)$se)
  #n.terms <- length(unlist(model.tables(x, type="effects", se=T)$se))
  report <- list(0)
  for(k in 2:n.terms){
  name.term <- names(model.tables(x, type="means")$n)[k-1]
  fitted  <- model.tables(x, type="means")$tables[[k]]
  num.means <- length(fitted)
  nis <- as.vector(model.tables(x, type="means")$n[[k-1]])
  sem <- se.terms[k-1]
  if(length(nis) == 1) nis <- rep(nis, num.means)
    #df <- x[[k]]$df.residual
    MSE <- deviance(x[[k]])/df 
    cat(paste("term", k, df, MSE, sqrt(MSE/df), se.terms[k-1], "\n", sep=" "))   
  #names(fitted)<-paste(c("M"),1:length(fitted),sep="")
  if(decreasing==TRUE){
  fitted<-sort(fitted, decreasing=TRUE)
  }
  else{
  fitted<-sort(fitted, decreasing=FALSE)
  }
   
  diffs<-matrix(0,num.means,num.means)
  for(i in 1:num.means){
    for(j in 1:num.means){
      diffs[i,j]<-fitted[j]-fitted[i]
    }
  }

  
 if (type=="LSD") {crit.diff<-LSD(nis, df, MSE, conf.level)}
  if (type=="NewmanKeuls"){crit.diff<-NewmanKeuls(nis, df, MSE, conf.level)} 
  if (type=="Duncan"){crit.diff<-Duncan(nis, df, MSE, conf.level)}
  if (type=="TukeyHSD"){crit.diff<-TukHSD(nis, df, MSE, conf.level)}
  if (type=="TukeyMRT"){crit.diff<-TukMRT(nis, df, MSE, conf.level)}
  if (type=="Scheffe"){crit.diff<-Scheffe(nis, df, MSE, conf.level)}
  if (type=="Bonferroni"){crit.diff<-Bonferroni(nis, df, MSE, conf.level)}
  if (type=="Sidak"){crit.diff<-Sidak(nis, df, MSE, conf.level)}
  
  sign<-matrix("",nrow=num.means,ncol=num.means)
  ext.lett<-c(letters,LETTERS)
  lettera<-1
  col<-1
  pos<-0
  posold<-0
  while(col<=num.means){
    
    if(abs(diffs[num.means,col])<=crit.diff[num.means,col]){
    
       for(i in col:num.means){
          sign[i,col]<-ext.lett[lettera]          
       }
      break                                                       
    }
    else{
      pos<-num.means
      while(abs(diffs[pos,col])>=crit.diff[pos,col]){        
        pos<-pos-1
        if (pos==col) break       
        }
      if (pos!=posold){
        for(i in col:pos){
          sign[i,col]<-ext.lett[lettera]
        }
        lettera<-lettera+1
        posold<-pos
      }
    
    }
   col<-col+1 
    
  }
  
  sign.lett<-sapply(1:num.means,function(nc)paste(sign[nc,],collapse=''))
  result<-data.frame(medie=fitted, test=sign.lett, row.names=names(fitted))  
#  cat(paste("term: ", name.term, "\n"))
#  print(result)
#  cat(paste("\n"))
#  print(paste("Critical.Difference = ", sqrt(2*MSE/nis), sep=" "))
  report[[k-1]] <- list(Effect=name.term, sem=sqrt(MSE/nis), sed=sqrt(2*MSE/nis), MCT=result, "Critical.Differences"=crit.diff[,1]) 
}
invisible(report)
}   
}
  
LSD <- function(nis, df, MSE, conf.level){
  T <- qt(conf.level/2, df, lower.tail=FALSE)
  min.diffs<-matrix(0,length(nis),length(nis))
  for(i in 1:length(nis)){
    for(j in 1:length(nis)){
    ifelse(df!=0, min.diffs[i,j]<-T*sqrt(MSE*(nis[i] + nis[j])/(nis[i]*nis[j])), min.diffs[i,j]<-T*MSE*sqrt(2))
    }
  }
  return(min.diffs)
}

NewmanKeuls <- function(nis, df, MSE, conf.level){
  min.diffs<-matrix(0,length(nis),length(nis))
  for(j in 1:(length(nis)-1)){
    for(i in (j+1):(length(nis))){
      pos<-i-j+1
      T <- qtukey(1-conf.level, pos, df)
      min.diffs[i,j]<-T*sqrt(MSE*(1/nis[i] + 1/nis[j])/2)
    }
  }
  min.diffs
}

Duncan <- function(nis, df, MSE, conf.level){
  min.diffs<-matrix(0,length(nis),length(nis))
  for(j in 1:(length(nis)-1)){
    for(i in (j+1):(length(nis))){
      pos<-i-j+1
      T <- qtukey(((1-conf.level)^(pos-1)), pos, df)
      min.diffs[i,j]<-T*sqrt(MSE*(1/nis[i] + 1/nis[j])/2)
    }
  }
  min.diffs
}

TukHSD <- function(nis, df, MSE, conf.level){
  min.diffs<-matrix(0,length(nis),length(nis))
  for(j in 1:(length(nis)-1)){
    for(i in (j+1):(length(nis))){
      T <- qtukey(1-conf.level, length(nis), df)
      min.diffs[i,j]<-T*sqrt(MSE*(1/nis[i] + 1/nis[j])/2)
    }
  }
  min.diffs
}

TukMRT <- function(nis, df, MSE, conf.level=.05){
  min.diffs<-matrix(0,length(nis),length(nis))
  for(j in 1:(length(nis)-1)){
    for(i in (j+1):(length(nis))){
      pos<-i-j+1
      T1 <- qtukey(1-conf.level, length(nis), df)
      T2 <- qtukey(1-conf.level, pos, df)
      T<-0.5*(T1+T2)
      min.diffs[i,j]<-T*sqrt(MSE*(1/nis[i] + 1/nis[j])/2)
    }
  }
  min.diffs
}

Scheffe <- function(nis, df, MSE, conf.level){
  r<-length(nis)
  min.diffs<-matrix(0,length(nis),length(nis))
  for(j in 1:(length(nis)-1)){
    for(i in (j+1):(length(nis))){
      T <- sqrt((r-1)*qf(1-conf.level,r-1,df))
      min.diffs[i,j]<-T*sqrt(MSE*(1/nis[i] + 1/nis[j]))
    }
  }
  min.diffs
}
Bonferroni <- function(nis, df, MSE, conf.level=.05){
  r<-length(nis)
  min.diffs<-matrix(0,length(nis),length(nis))
  for(j in 1:(length(nis)-1)){
    for(i in (j+1):(length(nis))){
      T <- qt(1-(conf.level)/(r*(r-1)),df)
      min.diffs[i,j]<-T*sqrt(MSE*(1/nis[i] + 1/nis[j]))
    }
  }
  min.diffs
}

Sidak <- function(nis, df, MSE, conf.level=.05){
  r<-length(nis)
  g <- r*(r-1)/2
  conf.level <- 1-(1-conf.level/2)^(1/g)
  min.diffs<-matrix(0,length(nis),length(nis))
  for(j in 1:(length(nis)-1)){
    for(i in (j+1):(length(nis))){
      T <- qt(1-(conf.level),df)
      min.diffs[i,j]<-T*sqrt(MSE*(1/nis[i] + 1/nis[j]))
    }
  }
  min.diffs
}
