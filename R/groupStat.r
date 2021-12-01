#crea le statistiche descrittive
groupStat <- function (x, cols, treat, FUN=sum){
    #x un data frame
    #cols il range di colonne da analizzare
    #treat il range di colonne factor
    means <- c()
    classe <- c()
    for(i in 1:length(treat)) classe<-paste(classe,x[,treat[i]], sep="_")
    for(i in 1:length(cols)) {
        stat <- c()
        if(is.numeric(x[,cols[i]])==TRUE) stat <- as.vector(tapply(x[,cols[i]], classe, FUN))
        means <- cbind(means, stat)
    }
    tesi<-c()
    for(i in 1:length(levels(factor(classe)))) {temp <- as.vector(strsplit(levels(factor(classe))[i], "_")[[1]])
    print(as.vector)
    tesi <- rbind(tesi, temp)}
   # print(means);print(dimnames(means)[[2]]); print(dimnames(x)[[2]])
 dimnames(means)[[2]] <- dimnames(x[cols])[[2]]
 row.names(means) <- levels(factor(classe))
 print(head(tesi))
 print(head(means))
 #means <- data.frame(tesi, means)
 #print(head(means))
 

 return(list(tesi, means))
}
