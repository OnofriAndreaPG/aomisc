plotnls <- function(fm, type = "means", xlim = NULL, res = 100, 
                    which = 3, ...){

    if (!inherits(fm, "nls"))
    stop("use only with \"nls\" objects")
    
    if(!is.numeric(which) || any(which < 1) || any(which > 3))
    stop("'which' must be in 1:3")
  
  if(which == 3){
    df <- eval(fm$data)
    mm <- fm$m
    cc <- fm$call
    pnms <- names(mm$getPars())
    form <- cc$formula
    rhsnms <- all.vars(form[[3]])
    vnms <- rhsnms[!(rhsnms %in% pnms)]
    if (length(vnms) > 1)
            stop("plotnls not yet implemented for >1 covariate")
    namList <- list(x = as.name(vnms), y = form[[2]])
    x <- df[,as.character(namList$x)]
    y <- df[,as.character(namList$y)]
    
    if(type == "means"){
                y <- tapply(y, list(factor(x)), mean)
                x <- tapply(x, list(factor(x)), mean)
                }
    
    if(is.null(xlim)) {
      xmin <- min(x)
      xmax <- max(x)  
    } else {
      xmin <- xlim[1]
      xmax <- xlim[2]
    } 
    xseq <- seq(xmin, xmax, 0.1)
    xseqDf <- data.frame(xseq)
    names(xseqDf) <- as.character(namList$x)
    newData <- predict(fm, newdata = xseqDf)
    
    plot(y ~ x, data = df, ...)
    points(newData ~ xseq, type = "l", col = "red")
  } else {
    #Get values
    x <- fm
    r <- residuals(x)
    yh <- predict(x) # != fitted() for glm
    w <- weights(x)
    
    if(!is.null(w)) { # drop obs with zero wt: PR#6640
      wind <- w != 0
      r <- r[wind]
      yh <- yh[wind]
      w <- w[wind]
      labels.id <- labels.id[wind]
    }
    
    n <- length(r)
    if(which == 2){
      s <- sqrt(deviance(x)/df.residual(x))
      ylab23 <- "Standardized residuals"
      r.w <- if (is.null(w)) r else sqrt(w) * r
      rs <- r.w/s
      }
    
    l.fit <- "Fitted values"
    
    ##---------- Do the individual plots : ----------
    if (which == 1) {
      ylim <- range(r, na.rm=TRUE)
      ylim <- extendrange(r = ylim, f = 0.08)
      plot(yh, r, xlab = l.fit, ylab = "Residuals", main = "Residuals vs Expected",
           ylim = ylim, type = "p", ...)
      abline(h = 0, lty = 3, col = "black")
    }
    if (which == 2) { ## Normal
      ylim <- range(rs, na.rm=TRUE)
      ylim[2L] <- ylim[2L] + diff(ylim) * 0.075
      # x.coord <- qnorm(ppoints(y))
      # y.coord <- scale(y, scale = T) 
      # plot(y.coord ~ x.coord, main = "Manual QQ-plot",
      # xlab = "Theoretical quantiles", ylab = "Standardised values")
      qq <- qqnorm(rs, main = "Normal Q-Q plot", ylab = ylab23, ylim = ylim, ...)
      qqline(rs, lty = 3, col = "gray50")
    }
  }  
    
}

