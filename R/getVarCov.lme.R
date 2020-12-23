getVarCov.lme <- function (obj, individuals, type = c("random.effects", "conditional", "marginal"), ...) {
    #Revised according to James Pustejovsky and AO (12/7/19)
    # ( see: https://www.jepusto.com/bug-in-nlme-getvarcov/)
    type <- match.arg(type)
    if (any("nlme" == class(obj))) 
        stop("not implemented for \"nlme\" objects")
    if (length(obj$group) > 1) 
        stop("not implemented for multiple levels of nesting")
    sigma <- obj$sigma
    D <- as.matrix(obj$modelStruct$reStruct[[1]]) * sigma^2
    
    if (type == "random.effects") {
        result <- D
    }
    else {
        result <- list()
        groups <- sort(obj$groups[[1]])
        ugroups <- unique(groups)
        if (missing(individuals)) 
            individuals <- as.matrix(ugroups)[1, ]
        if (is.numeric(individuals)) 
            individuals <- ugroups[individuals]
        for (individ in individuals) {
            indx <- which(individ == ugroups)
            if (!length(indx)) 
                stop(gettextf("individual %s was not used in the fit", 
                  sQuote(individ)), domain = NA)
            if (is.na(indx)) 
                stop(gettextf("individual %s was not used in the fit", 
                  sQuote(individ)), domain = NA)
            ind <- groups == individ
            if (!is.null(obj$modelStruct$corStruct)) {
                V <- corMatrix(obj$modelStruct$corStruct)[[as.character(individ)]]
            }
            else V <- diag(sum(ind))
            if (!is.null(obj$modelStruct$varStruct)) 
                sds <- 1/varWeights(obj$modelStruct$varStruct)[ind]
            else sds <- rep(1, sum(ind))
            sds <- obj$sigma * sds
            cond.var <- t(V * sds) * sds
            dimnames(cond.var) <- list(1:nrow(cond.var), 1:ncol(cond.var))
            if (type == "conditional") 
                result[[as.character(individ)]] <- cond.var
            else {
                Z <- model.matrix(obj$modelStruct$reStruc, getData(obj))[ind, 
                  , drop = FALSE]
                result[[as.character(individ)]] <- cond.var + 
                  Z %*% D %*% t(Z)
            }
        }
    }
    class(result) <- c(type, "VarCov")
    attr(result, "group.levels") <- names(obj$groups)
    result
}
