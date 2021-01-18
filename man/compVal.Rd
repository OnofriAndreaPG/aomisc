\name{compVal}
\alias{compVal}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Pairwise comparisons of 
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
compVal(object, SE, parameterNames, operator = "-", df.residual = NA, method = "pairwise", adjust = "holm", decreasing = FALSE, level = 0.05)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{
%%     ~~Describe \code{object} here~~
}
  \item{SE}{
%%     ~~Describe \code{SE} here~~
}
  \item{parameterNames}{
%%     ~~Describe \code{parameterNames} here~~
}
  \item{operator}{
%%     ~~Describe \code{operator} here~~
}
  \item{df.residual}{
%%     ~~Describe \code{df.residual} here~~
}
  \item{method}{
%%     ~~Describe \code{method} here~~
}
  \item{adjust}{
%%     ~~Describe \code{adjust} here~~
}
  \item{decreasing}{
%%     ~~Describe \code{decreasing} here~~
}
  \item{level}{
%%     ~~Describe \code{level} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (object, SE, parameterNames, operator = "-", df.residual = NA, 
    method = "pairwise", adjust = "holm", decreasing = FALSE, 
    level = 0.05) 
{
    strVal <- factor(parameterNames, levels = as.character(parameterNames))
    presentVec <- strVal
    lenPV <- length(presentVec)
    if (lenPV < 2) {
        stop("No parameters to compare")
    }
    parm <- as.vector(object)
    varMat <- diag(SE^2)
    if (identical(operator, "/")) {
        hypVal <- 1
        fct <- function(ind) {
            parm[ind[1]]/parm[ind[2]]
        }
        dfct <- function(ind) {
            transVec <- c(1/parm[ind[2]], -parm[ind[1]]/(parm[ind[2]]^2))
            sqrt(transVec \%*\% varMat[ind, ind] \%*\% transVec)
        }
    }
    if (identical(operator, "-")) {
        hypVal <- 0
        fct <- function(ind) {
            parm[ind[1]] - parm[ind[2]]
        }
        transVec <- c(1, -1)
        dfct <- function(ind) {
            sqrt(transVec \%*\% varMat[ind, ind] \%*\% transVec)
        }
    }
    lenRV <- lenPV * (lenPV - 1)/2
    cpMat <- matrix(0, lenRV, 4)
    compParm <- rep("", lenRV)
    degfree <- ifelse(is.na(df.residual) == T, Inf, df.residual)
    pFct <- function(x) {
        pt(x, degfree)
    }
    strParm <- strVal
    k <- 1
    for (i in 1:lenPV) {
        for (j in 1:lenPV) {
            if (j <= i) {
                next
            }
            cpMat[k, 1] <- fct(presentVec[c(i, j)])
            cpMat[k, 2] <- dfct(presentVec[c(i, j)])
            tVal <- (cpMat[k, 1] - hypVal)/cpMat[k, 2]
            cpMat[k, 3] <- tVal
            cpMat[k, 4] <- pFct(-abs(tVal)) + (1 - pFct(abs(tVal)))
            compParm[k] <- paste(strParm[presentVec[c(i, j)]], 
                collapse = operator)
            k <- k + 1
        }
    }
    cpMat <- data.frame(cpMat)
    row.names(cpMat) <- compParm
    colnames(cpMat) <- c("Estimate", "Std. Error", "t-value", 
        "p-value")
    adjusted.P <- p.adjust(cpMat$`p-value`, method = as.character(adjust))
    cpMat$`p-value` <- as.vector(adjusted.P)
    if (method == "pairwise") {
        return(cpMat)
    }
    else {
        p.logic <- as.vector(adjusted.P)
        p.logic <- ifelse(p.logic > 0.05, FALSE, TRUE)
        names(p.logic) <- sub("/", "-", as.vector(row.names(cpMat)))
        tmp <- data.frame(parNam = parameterNames, obj = object)
        if (decreasing == T) {
            Letters <- multcompView::multcompLetters3("parNam", 
                "obj", p.logic, data = tmp, reverse = F, threshold = level)
            parMat <- data.frame(Value = object, SE = SE)
            row.names(parMat) <- parameterNames
            parMat <- parMat[order(-parMat$Value), ]
        }
        else {
            Letters <- multcompView::multcompLetters3("parNam", 
                "obj", p.logic, data = tmp, reverse = T, threshold = level)
            parMat <- data.frame(Value = object, SE = SE)
            row.names(parMat) <- parameterNames
            parMat <- parMat[order(-parMat$Value), ]
        }
        parMat$CLD = as.character(Letters$Letters)
        if (decreasing == F) {
            parMat <- parMat[order(parMat$Value), ]
        }
        else {
            parMat <- parMat[order(-parMat$Value), ]
        }
        return(parMat)
    }
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory (show via RShowDoc("KEYWORDS")):
% \keyword{ ~kwd1 }
% \keyword{ ~kwd2 }
% Use only one keyword per line.
% For non-standard keywords, use \concept instead of \keyword:
% \concept{ ~cpt1 }
% \concept{ ~cpt2 }
% Use only one concept per line.
