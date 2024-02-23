biplot.AMMIobject <- function(x, biplot = 1, xlim=NULL, ylim=NULL, 
                              elabels=NULL, glabels=NULL, quad=FALSE, 
                              cexG=0.9, cexE=0.9, xlab=NULL, ylab=NULL, 
                              font=1, ...){
    # Edited on 26/5/2023 
    obj <- x
    E <- obj$environment_scores
    G <- obj$genotype_scores
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
            labels = if(is.null(elabels)) sub(' +$', '', as.vector(names(envir.mean))) else elabels, 
            #adj = c(0.5, 0.5), 
            col = "red", font=font, 
            cex = cexE)
        points(var.mean, G[, 1], "n", col = "blue", lwd = 5)
        text(var.mean, G[, 1], 
            labels = if(is.null(glabels)) sub(' +$', '', as.vector(names(var.mean))) else glabels, 
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
        if (quad == TRUE){
        abline(h = 0, lty = 5)
        abline(v = 0, lty=5)
        }
        # if (percDev != FALSE){
        #   expvar <- c(paste(round(obj$mult_Interaction$Perc_of_Total_SS[1]),"%"), paste(round(obj$mult_Interaction$Perc_of_Total_SS[2]),"%"))
        #   arrows(percDev[1], percDev[2], (percDev[1]+max(c(E[, 1], G[, 1]))/7), percDev[2], length=0.15)
        #   text((percDev[1]+max(c(E[, 1], G[, 1]))/7), percDev[2], label=expvar[1], pos=4)
        #   arrows(percDev[1], percDev[2],percDev[1], (percDev[2]+max(c(E[, 1], G[, 1]))/7), length=0.15)
        #   text(percDev[1], (percDev[2]+max(c(E[, 1], G[, 1]))/7), label=expvar[2], pos=3)
        # }
    }
}

biplot.GGEobject <- function(x, biplot = 1, xlim=NULL, ylim=NULL, 
                              elabels=NULL, glabels=NULL, quad=FALSE, 
                              cexG=0.9, cexE=0.9, xlab=NULL, ylab=NULL, 
                              font=1, ...){
  obj <- x
  biplot <- 2
  biplot.AMMIobject(obj, biplot, xlim, ylim, elabels, glabels, 
                    quad, cexG, cexE, xlab, ylab, font)
  
}



