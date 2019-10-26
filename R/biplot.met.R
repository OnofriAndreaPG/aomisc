biplot.met <- function(obj, biplot = 1, xlim=NULL, ylim=NULL, elabels=NULL, 
                       glabels=NULL, quad=F, percDev=F, size=2, cexG=0.9, cexE=0.9,
                       xlab=NULL, ylab=NULL, font=1){
    if (class(obj) == "GGEobject") biplot <- 2                  
    E <- obj$environment_scores
    G <- obj$genotype_scores
    print(G)
    envir.mean <- obj$environment_means
    var.mean <- obj$genotype_means
    overall.mean <- mean(envir.mean)
    if (biplot == 1) {
        plot(1, type = "n", 
            xlim = if(is.null(xlim)) range(c(envir.mean, var.mean)) else xlim, 
            ylim = if(is.null(ylim)) range(c(E[, 1], G[, 1])) else ylim, 
            xlab = if(is.null(xlab)) "VarY" else xlab, 
            ylab = if(is.null(ylab)) "PC 1" else ylab,)
        points(envir.mean, E[, 1], "n", col = "red", lwd = 5)
        text(envir.mean, E[, 1],
            labels = if(is.null(elabels)) sub(' +$', '', as.vector(row.names(envir.mean))) else elabels, 
            #adj = c(0.5, 0.5), 
            col = "red", font=font, 
            cex = cexE)
        points(var.mean, G[, 1], "n", col = "blue", lwd = 5)
        text(var.mean, G[, 1], 
            labels = if(is.null(glabels)) sub(' +$', '', as.vector(row.names(var.mean))) else glabels, 
            #adj = c(0.5, 0.5), 
            col = "blue", font=font, 
            cex = cexG)
        abline(h = 0, lty = 5)
        abline(v = overall.mean, lty=5) 
    }
    else {
        plot(1, type = "n", 
        xlim = if(is.null(xlim)) range(c(E[, 1], G[, 1])) else xlim, 
        ylim = if(is.null(ylim)) range(c(E[, 2], G[, 2])) else ylim,
        xlab = "PC 1", ylab = "PC 2")
        points(E[, 1], E[, 2], "n", col = "red", lwd = 5)
        text(E[, 1], E[, 2], 
        	labels = if(is.null(elabels)) sub(' +$', '', as.vector(row.names(E))) else elabels,
        	col = "red", font=2, 
        	cex=cexE)
        points(G[, 1], G[, 2], "n", col = "blue", lwd = 5)
        text(G[, 1], G[, 2], 
        labels = if(is.null(glabels)) sub(' +$', '', as.vector(row.names(G))) else glabels, 
        	#adj = c(0.5, 0.5), 
        	col = "blue", font=font, 
        	cex=cexG)
        if (quad==TRUE){
        abline(h = 0, lty = 5)
        abline(v = 0, lty=5)
        }
        if (percDev != FALSE){
          expvar <- c(paste(round(obj$GGE_summary$Perc_of_Total_SS[1]),"%"), paste(round(obj$GGE_summary$Perc_of_Total_SS[2]),"%"))
          arrows(percDev[1], percDev[2], (percDev[1]+max(c(E[, 1], G[, 1]))/7), percDev[2], length=0.15)
          text((percDev[1]+max(c(E[, 1], G[, 1]))/7), percDev[2], label=expvar[1], pos=4)
          arrows(percDev[1], percDev[2],percDev[1], (percDev[2]+max(c(E[, 1], G[, 1]))/7), length=0.15)
          text(percDev[1], (percDev[2]+max(c(E[, 1], G[, 1]))/7), label=expvar[2], pos=3)
        }
    }
}

segment <- function(x1, y1, x2, y2, lty=NULL){
  curve(((x - x1)*(y2 - y1)) / (x2 - x1) + y1, from=x1, to=x2, lty = if(is.null(lty)) 1 else lty, add=T)
}

perp2p <- function(x1, y1, x2, y2, x3, y3, lty=NULL){
  m <- (y2 - y1)/(x2 - x1)
  curve(-1/m*(x-x3)+y3, from=if(max(x1,x2)>=x3 & min(y1,y2)<=y3) x3, to=if(min(x1,x2)<=x3 & max(y1,y2)>=y3) x3, lty = if(is.null(lty)) 1 else lty, add=T) 
}
biplot.polygon <- function(vertex, obj){
E <- obj$Environment_scores
G <- obj$Genotype_scores
  for (i in 1:(length(vertex)-1)){
  segment(G[vertex[i],1], G[vertex[i],2], G[vertex[i+1],1], G[vertex[i+1],2], lty=2)
  #perp2p(G[(vertex[i]),1], G[(vertex[i]),2], G[(vertex[i+1]),1], G[(vertex[i+1]),2], 0, 0)
}
  segment(G[vertex[1],1], G[vertex[1],2], G[vertex[length(vertex)],1], G[vertex[length(vertex)],2], lty=2)
  #perp2p(G[vertex[1],1], G[vertex[1],2], G[vertex[length(vertex)],1], G[vertex[length(vertex)],2], 0, 0)
}
