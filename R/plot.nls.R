plot.nls <-
  function (x, which = 1, 
            caption = "Residuals vs Fitted",
            sub.caption = NULL, main = "", ...,
            label.pos = c(4,2), cex.caption = 1, cex.oma.main = 1.25)
  {
    if (!inherits(x, "nls"))
      stop("use only with \"nls\" objects")
    if(!is.numeric(which) || any(which < 1) || any(which > 2))
      stop("'which' must be in 1:2")
    
    
    #Get values
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
      qq <- qqnorm(rs, main = "Normal Q-Q plot", ylab = ylab23, ylim = ylim, ...)
      qqline(rs, lty = 3, col = "gray50")
    }
  }