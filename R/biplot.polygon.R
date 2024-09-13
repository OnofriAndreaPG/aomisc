biplot.polygon <- function(obj, vertex){
  # It draws a which-won-where pattern on a GGE biplot
  # it should return an error if there is no existing plot
  if(dev.cur() == 1L ) stop("There is no existing plot")
  if(!inherits(obj, "GGEobject")) stop("The polygon must be superimposed to a GGE biplot")
  
  segment <- function(x1, y1, x2, y2, lty=NULL, col = 1){
  # This function draws a segment joining two points.
  curve(((x - x1)*(y2 - y1)) / (x2 - x1) + y1, 
        from=x1, to=x2, lty = if(is.null(lty)) 1 else lty, add=T, col = col)
  }
  perp2p <- function(x1, y1, x2, y2, x3, y3, lty=NULL, col = 1){
  # It draws the perpendicular to a segment from P1(x1,y1) to P2(x2, y2)
  # passing through the point P3(x3,y3)
  m <- (y2 - y1)/(x2 - x1)
  # fromP <- if(max(x1,x2)>=x3 & min(y1,y2)<=y3) x3 else NULL
  if(all(c(x1,x2) >=x3)){
    fromP <- x3
    toP <- NULL
  } else if(all(c(x1,x2) <= x3)){
    fromP <- NULL
    toP <- x3 
  } else if (all(c(y1,y2) >= y3)){
    fromP <- NULL
    toP <- x3
  } else {
    fromP <- NULL
    toP <- x3
  }
  
  # toP <- if(min(x1,x2)<=x3 & max(y1,y2)>=y3) x3 else NULL
  # toP <- if(all(c(x1,x2)<=x3) ) x3 else NULL
  
  curve(-1/m*(x - x3) + y3, from = fromP, to=toP, 
        lty = if(is.null(lty)) 1 else lty, add=T, col = col)
  }

  E <- obj$environment_scores
  G <- obj$genotype_scores
    for (i in 1:(length(vertex)-1)){
    segment(G[vertex[i],1], G[vertex[i],2], G[vertex[i+1],1], G[vertex[i+1],2], lty=2)
    perp2p(G[(vertex[i]),1], G[(vertex[i]),2], G[(vertex[i+1]),1], G[(vertex[i+1]),2], 0, 0)
    }
    segment(G[vertex[1],1], G[vertex[1],2], G[vertex[length(vertex)],1], G[vertex[length(vertex)],2], lty=2)
    perp2p(G[vertex[1],1], G[vertex[1],2], G[vertex[length(vertex)],1], G[vertex[length(vertex)],2], 0, 0)
}