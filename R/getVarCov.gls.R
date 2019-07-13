getVarCov.gls <- function (obj, individual = 1, ...) {
    #Revised according to James Pustejovsky and AO (12/7/19)
    # ( see: https://www.jepusto.com/bug-in-nlme-getvarcov/)
    S <- corMatrix(obj$modelStruct$corStruct)[[individual]]
    if (!is.null(obj$modelStruct$varStruct)) {
        ind <- sort(obj$groups) == obj$groups[[individual]]
        vw <- 1 / varWeights(obj$modelStruct$varStruct)[ind]
    }
    else vw <- rep(1, nrow(S))
    vars <- (obj$sigma * vw)^2
    result <- t(S * sqrt(vars)) * sqrt(vars)
    class(result) <- c("marginal", "VarCov")
    attr(result, "group.levels") <- names(obj$groups)
    result
}
