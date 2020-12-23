# Serve per poter fare confronti passando un semplice
# vettore di medie ed errori standard
tukeyMat <- function(obj, lev = NULL) {
  # obj: vector of means, possibly named
  # lev: optional names for means
  # Crea una matrice di contrasti Tukey-type
    lfm <- diag(length(obj))
    nlev <- nrow(lfm)
    rn <- names(obj)
    a <- attr(lfm, "grid")
    if (is.null(lev)) {
        if (!is.null(a)) {
            lev <- apply(a, 1, paste, collapse = ":")
        } else if (!is.null(rn)) {
            lev <- rn
        } else {
            lev <- as.character(1:nlev)
        }
    }
    cbn <- utils::combn(seq_along(lev), 2)
    M <- lfm[cbn[1, ], ] - lfm[cbn[2, ], ]
    if (is.vector(M)) {
        dim(M) <- c(1, length(M))
    }
    rownames(M) <- paste(lev[cbn[1, ]], lev[cbn[2, ]], sep = "-")
    return(M)
}

# Serve per ottenere le lettere in modo semplice
# Passando solo un glht object
cld2 <- function(obj, pval = 0.05){
  if(class(obj) != "glht"){
    print("this function is only implemented for glht objects")
    stop()
  }
  gh <- summary(obj)
  p.logic <- as.vector(gh$test$pvalues)
  p.logic <- ifelse(p.logic > pval, FALSE, TRUE)
  names(p.logic) <- sub(" - ", "-", as.vector(dimnames(gh$linfct)[[1]]))
  Letters <- multcompView::multcompLetters(p.logic)
  return(Letters)
}

pairComp <- function(parm, SE, nams){
  names(parm) <- nams
  pairList <- list(coef = parm, vcov = diag(SE^2), df = Inf)
  class(pairList) = "parm"
  gh <- glht(pairList, linfct = tukeyMat(parm))
  lett <- cld2(gh)
  letters <- data.frame(Mean = parm, SE = SE, CLD = lett$Letters)
  letters$Zval = letters$Mean/letters$SE
  letters$prob = 1 - pnorm(letters$Zval)
  returnList = list(pairs = summary(gh), Letters = letters)
  return(returnList)
}

contrMat <- function(n, type = c("Dunnett", "Tukey", "Sequen", "AVE",
                                 "Changepoint", "Williams", "Marcus",
                                 "McDermott", "UmbrellaWilliams", "GrandMean",
                                 "Mean"), 
                     base = 1) {

    if (length(n) < 2) stop("less than two groups")
    if (!is.numeric(n)) stop(sQuote("n"), " is not numeric")
    m <- NULL
    type <- match.arg(type)
    if (type %in% c("AVE", "Williams", "McDermott") && length(n) < 3)
        stop("less than three groups")
    k <- length(n)
    if (base < 1 || base > k) stop("base is not between 1 and ", k)
    CM <- c()
    rnames <- c()
    if (!is.null(names(n)))
        varnames <- names(n)
    else 
        varnames <- 1:length(n)

    kindx <- 1:k

    switch(type, 
      "Dunnett" = {
        for(i in kindx[-base])
            CM <- rbind(CM, as.numeric(kindx == i) - as.numeric( kindx == base))
        rnames <- paste(varnames[kindx[-base]], "-", varnames[base])
    }, "Tukey" = {
        for (i in 1:(k-1)) {
            for(j in (i+1):k) {
                CM  <- rbind(CM, as.numeric(kindx==j)-as.numeric(kindx==i))
                rnames <- c(rnames, paste(varnames[j], "-", varnames[i]))
            }
        }
    }, "Sequen" =  {
        for (i in 2:k) {
            CM  <- rbind(CM, as.numeric(kindx==i)-as.numeric(kindx==i-1))
            rnames <- c(rnames, paste(varnames[i], "-", varnames[i-1]))
        }
    }, "AVE" = {
        help <- c(1,  -n[2:k]/sum(n[2:k]))
        CM <- rbind(CM, help)
        for (i in 2:(k-1)) {
            x <- sum(n[1:(i-1)])+sum(n[(i+1):k])
            help <- c(-n[1:(i-1)]/x, 1, -n[(i+1):k]/x)
            CM <- rbind(CM, help)
        }
        help <- c(-n[1:(k-1)]/sum(n[1:(k-1)]), 1)
        CM  <- rbind(CM, help)
        rnames <- paste("C", 1:nrow(CM))
    }, "Changepoint" = {
        for (i in 1:(k-1)) {
            help <- c(-n[1:i]/sum(n[1:i]), n[(i+1):k]/sum(n[(i+1):k]))
            CM <- rbind(CM, help)
        }
        rnames <- c(rnames, paste("C", 1:nrow(CM)))
    }, "Williams" = {
        for (i in 1:(k-2)) {
            help <-  c(-1, rep(0, k-i-1), n[(k-i+1):k]/sum(n[(k-i+1):k]))
            CM <- rbind(CM, help)
        }
        help <- c(-1, n[2:k]/sum(n[2:k]))
        CM <- rbind(CM, help)
        rnames <- c(rnames, paste("C", 1:nrow(CM)))
    }, "Marcus" = {
        cm1 <- matrix(0, nrow=k-1, ncol=k)
        cm2 <- cm1
        for (i in 1:(k-1)) {
            cm1[i,(i+1):k] <- n[(i+1):k]/sum(n[(i+1):k])
            cm2[i,1:i] <- n[1:i]/sum(n[1:i])
        }
        ### row <- k*(k-1)/2
        index <- 1
        for (i in 1:(k-1)) {
            for (j in 1:i) {
                help <- cm1[i,]-cm2[j,]
                CM <- rbind(CM, help)
                index <- index+1
            }
        }
        rnames <- c(rnames, paste("C", 1:nrow(CM)))
     }, "McDermott" = {
         for(i in 1:(k-2)) {
             help  <- c(-n[1:i]/sum(n[1:i]), 1, rep(0, k-i-1))
             CM <- rbind(CM, help)
         }
         help <- c(-n[1:(k-1)]/sum(n[1:(k-1)]), 1)
         CM  <- rbind(CM, help)
         rnames <- c(rnames, paste("C", 1:nrow(CM)))
    }, "Tetrade" = {
        if (is.null(m)) stop(sQuote("m"), " is missing")
        a <- length(n)
        b <- length(m)
        if (!is.null(names(m)))
            varnamesm <- names(m)
        else 
            varnamesm <- 1:length(m)
	idi <- 1:a
	idj <- 1:b
        for (i1 in 1:(a-1)) {
            for (i2 in (i1+1):a) {
	        for (j1 in 1:(b-1)) {
        	    for (j2 in (j1+1):b) {
                	CM <- rbind(CM, kronecker( ( as.numeric(idi==i1)-as.numeric(idi==i2) ),
                                                   ( as.numeric(idj==j1)-as.numeric(idj==j2) ) ) ) 
		        rnames <- c(rnames, paste( "(", paste(varnames[i1], varnamesm[j1], sep = ":"), "-", 
                                                        paste(varnames[i1], varnamesm[j2], sep = ":"), ")", "-", 
                                                   "(", paste(varnames[i2], varnamesm[j1], sep = ":"), "-", 
                                                        paste(varnames[i2], varnamesm[j2], sep = ":"), ")",  sep=""))
            	    }
        	}
	    }
        }
    }, "UmbrellaWilliams" = {
        for (j in 1:(k-1)) {
            for (i in 1:(k - j)) {
                helper <- c(-1, rep(0, k - i - j),
                    n[((k - i + 1):k)-(j-1)]/sum(n[((k - i + 1):k)-(j-1)]),
                    rep(0, j-1))
                CM <- rbind(CM, helper)
            }
        }
        rnames <- c(rnames, paste("C", 1:nrow(CM)))
    }, "GrandMean" = {
        CM <- matrix(rep(-n/sum(n), k), nrow = k, byrow = TRUE)
        diag(CM) <- diag(CM) + 1
        rnames <- c(rnames, paste("C", 1:nrow(CM)))
    }, "Mean" = {
        CM <- matrix(0, nrow = k, ncol = k, byrow = TRUE)
        diag(CM) <- 1
        rnames <- c(rnames, paste("C", 1:nrow(CM)))
    })

    rownames(CM) <- rnames
    if (type == "Tetrade")
      colnames(CM) <- NULL ###levels(interaction(varnames, varnamesm))
    else 
      colnames(CM) <- varnames
    attr(CM, "type") <- type
    class(CM) <- c("contrMat", "matrix")
    CM
}

mcp2matrix <- function(model, linfct) {

    ### extract factors and contrasts
    fc <- factor_contrasts(model)
    contrasts <- fc$contrasts
    factors <- fc$factors
    intercept <- fc$intercept
    mf <- fc$mf
    mm <- fc$mm

    alternative <- NULL

    ### linear hypotheses
    if (!is.list(linfct) || is.null(names(linfct)))
        stop(sQuote("linfct"), "is not a named list")
    nhypo <- names(linfct)
    checknm <- nhypo %in% rownames(factors)
    if (!all(checknm)) 
        stop("Variable(s) ", sQuote(nhypo[!checknm]), " have been specified in ",
             sQuote("linfct"), " but cannot be found in ", sQuote("model"), "! ")
    if (any(checknm)) {
        checknm <- sapply(mf[nhypo[checknm]], is.factor)
        if (!all(checknm))
            stop("Variable(s) ", sQuote(paste(nhypo[!checknm], collapse = ", ")), " of class ", 
                  sQuote(paste(sapply(mf[nhypo[!checknm]], class), collapse = ", ")), 
                  " is/are not contained as a factor in ", sQuote("model"), ".")
    }
    m <- c()
    ctype <- c()
    for (nm in nhypo) {
        if (is.character(linfct[[nm]])) {

            Kchr <- function(kch) {
                ### check if kch is suitable as `type' argument to `contrMat'
                types <- eval(formals(contrMat)$type)
                pm <- pmatch(kch, types)
                ### if yes, compute K from `contrMat'
                if (!is.na(pm)) {
                    tmpK <- contrMat(table(mf[[nm]]), type = types[pm])
                    ctype <<- c(ctype, types[pm])
                } else {
                    ### if not, interpret kch as an expression
                    tmp <-  chrlinfct2matrix(kch, levels(mf[[nm]]))
                    tmpK <- tmp$K
                    m <<- c(m, tmp$m)
                    if (is.null(alternative)) {
                        alternative <<- tmp$alternative
                    } else {
                        if (tmp$alternative != alternative)
                            stop("mix of alternatives currently not implemented")
                    }
                }
                print(row.names(kch))
                if (is.null(rownames(tmpK)))
                    rownames(tmpK) <- paste(kch, 1:nrow(tmpK), sep = "_")
                if (length(nhypo) > 1)
                    rownames(tmpK) <- paste(nm, rownames(tmpK), sep = ": ")
                list(K = tmpK)
            }
            
            tmp <- lapply(linfct[[nm]], Kchr)
            linfct[[nm]] <- do.call("rbind", lapply(tmp, function(x) x$K))
        }
    }

    ### transform linear hypotheses using model contrasts
    hypo <- vector(mode = "list", length = length(nhypo))
    names(hypo) <- nhypo

    for (nm in nhypo) {
        ### extract contrast matrix for each factor from model fit
        if (is.character(contrasts[[nm]])) {
            C <- do.call(contrasts[[nm]], 
                         list(n = nlevels(mf[[nm]])))
        } else {
            C <- contrasts[[nm]]
        }
        ### and transform the original linear hypotheses 
        ### K beta to K C beta^* 
        if (intercept || (!intercept && nm != colnames(factors)[1])) {
            Kstar <- linfct[[nm]] %*% C
        } else {
            ### model.matrix has `contrasts' argument even if no intercept
            ### was fitted and the contrast actually hasn't been applied
            ### This is, however, only the case for the _first_ factor
            Kstar <- linfct[[nm]]
        }
        pos <- factors[nm,] == 1
        ### interaction terms (if any)
        if (sum(pos) > 1)
            warning("covariate interactions found -- ", 
                    "default contrast might be inappropriate")
        hypo[[nm]] <- list(K = Kstar,
            where = attr(mm, "assign") %in% which(nm == colnames(factors)))
    }

    ### combine all single matrices computed so far into
    ### one matrix of all linear hypoheses
    Ktotal <- matrix(0, nrow = sum(sapply(hypo, function(x) nrow(x$K))),
                     ncol = ncol(mm))
    colnames(Ktotal) <- colnames(mm)

    count <- 1
    for (h in hypo) {
        Ktotal[count:(count + nrow(h$K) - 1), h$where] <- h$K
        count <- count + nrow(h$K)
    }
    if (!is.matrix(Ktotal)) Ktotal <- matrix(Ktotal, nrow = 1)
    rownames(Ktotal) <- unlist(lapply(hypo, function(x) rownames(x$K)))

    if (is.null(ctype))
        ctype <- "User-defined"
    ctype <- paste(unique(ctype), collapse = ", ")
    attr(Ktotal, "type") <- ctype

    if (length(m) == 0) m <- 0
    list(K = Ktotal, m = m, alternative = alternative, type = ctype)
}

factor_contrasts <- function(model) {

    ### extract model matrix, frame and terms
    mm <- try(model.matrix(model))
    if (inherits(mm, "try-error"))
        stop("no ", sQuote("model.matrix"), " method for ", 
             sQuote("model"), " found!")

    mf <- try(model.frame(model))
    if (inherits(mf, "try-error"))
        stop("no ", sQuote("model.frame"), " method for ", 
             sQuote("model"), " found!")

    tm <- try(terms(model))
    if (inherits(tm, "try-error"))
        stop("no ", sQuote("terms"), " method for ", 
             sQuote("model"), " found!")

    list(contrasts = attr(mm, "contrasts"),
         factors = attr(tm, "factors"),
         intercept = attr(tm, "intercept") != 0,
         mm = mm, 
         mf = mf)
}

