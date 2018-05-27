makeSurv <- function(dati, tratt, numSemi, tempi) {
#Preparazione file. 
#dati = dataframe in cui, per ogni capsula petri (riga) e per ogni
#tempo (colonna) sono listati il numero dei semi germinati. I nomi di colonna
#sono i tempi di osservazione
#
#tratt = dataframe in cui, per ogni capsula Petri (riga) sono listati i livelli
#di ogni tesi sperimentale (colonna)
#
#numSemi = vettore che per ogni capsula Petri lista il numero di semi vitali,
#in ogni capsula Petri (num. iniziale - num. finale nonviable)
#
#tempi = vettore dei tempi di osservazione 

#Sostituisce i NAs con 0 e prepara il file
dati[is.na(dati)] <- 0

#crea i vettori giorni, cens, group;
giorni <- c(); cens <- c(); group <- data.frame()
final.date <- max(tempi)
for (j in 1:length(dati[,1])){ #Per ogni riga
    for(i in 1:length(dati[1,])) { #Per ogni colonna
        giorni <- c(giorni, rep(tempi[i], dati[j,i]))
        cens <- c(cens, rep(1, dati[j,i]))
    }
    #censorati
    giorni <- c(giorni, rep(final.date, (numSemi[j]-sum(dati[j,]))))
    cens <- c(cens, rep(0, (numSemi[j]-sum(dati[j,]))))
    #Nomi dei trattamenti
	 for(z in 1:numSemi[j]) group <- rbind(group, tratt[j,])
}
dimnames(group)[[2]] <- dimnames(tratt)[[2]]
rep(dimnames(dati)[[2]], each=numSemi[j])
datiFin <- data.frame(time=giorni, status=cens)
datiFin <- cbind(group, datiFin)
return(datiFin)

}

reorderFactor <- function(X, ordineRev) {
tratt2 <- X
for (i in 1:length(ordineRev)) tratt2 <- relevel(tratt2, ref=levels(tratt)[ordineRev[i]])
return(tratt2)
}