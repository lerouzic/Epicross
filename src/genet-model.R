
# Sets the coefficient loadings (additive, dominance, 
# add x add, add x dom, and dom x dom)
# for parental, F1, and F2 populations. 

model1 <- rbind(
	P1=c(a1=1,   a2=0,   d=0,   aa=0,    ad=0, dd=0),
	P2=c(a1=0,   a2=1,   d=0,   aa=0,    ad=0, dd=0), 
	F1=c(a1=1/2, a2=1/2, d=1,   aa=-1/2, ad=0, dd=1/2),
	F2=c(a1=1/2, a2=1/2, d=1/2, aa=-1/2, ad=0, dd=1/8))

model2 <- rbind(
	P=c(a=1,    d=0,   aa=0,    ad=0, dd=0),
	F1=c(a=1/2, d=1,   aa=-1/2, ad=0, dd=1/2),
	F2=c(a=1/2, d=1/2, aa=-1/2, ad=0, dd=1/8))
	
MODEL <- model2


# Design matrices
Z.2pop.F2 <- function(P1, P2, F1, F2) {
	nP1 <- length(P1)
	nP2 <- length(P2)
	nF1 <- length(F1)
	nF2 <- length(F2)
	list(
		mu = c(rep( 1, nP1), rep( 1, nP2), rep(1, nF1), rep(1, nF2)), 
		a  = c(rep( 1, nP1), rep(-1,nP2),  rep(0, nF1), rep(0, nF2)), 
		d  = c(rep(-1, nP1), rep(-1, nP2), rep(1, nF1), rep(0, nF2)), 
		aa = c(rep( 1, nP1), rep( 1, nP2), rep(0, nF1), rep(0, nF2)), 
		ad = c(rep(-1, nP1), rep( 1, nP2), rep(0, nF1), rep(0, nF2)), 
		dd = c(rep( 1, nP1), rep( 1, nP2), rep(0, nF1), rep(0, nF2)))
}

Z.Npop.Fmu <- function(fP1, fP2, fgen, grandmean=FALSE) {
	stopifnot(length(fP1) == length(fP2), length(fP1) == length(fgen))
	
	n <- length(fP1)                # Number of measurements
	f <- factor(unique(c(fP1,fP2))) # Populations
	fi <- paste(ifelse(fP1 < fP2, fP1, fP2), ifelse(fP1 < fP2, fP2, fP1), sep="x")
	i <- factor(unique(fi[fP1 != fP2]))
	
	mu <- rep(1, length(fP1))
	a  <- matrix(0, nrow=n, ncol=length(f), dimnames=list(c(), as.character(f)))
	d  <- matrix(0, nrow=n, ncol=length(i), dimnames=list(c(), as.character(i)))
	aa <- matrix(0, nrow=n, ncol=length(i), dimnames=list(c(), as.character(i)))
	dd <- matrix(0, nrow=n, ncol=length(i), dimnames=list(c(), as.character(i)))
	
	for (j in seq_along(fP1)) {
		fi1 <- which(f == fP1[j])
		fi2 <- which(f == fP2[j])
		
		a[j, fi1] <- a[j, fi1] + 1/2
		a[j, fi2] <- a[j, fi2] + 1/2
		
		if (fP1[j] != fP2[j]) {
			ii <- which(i == fi[j])
			aa[j,ii] <- -1
			if (fgen[j] == "F1") {
				d[j,ii]  <- 2
			} else {
				d[j,ii]  <- 1
				dd[j,ii] <- -1
			}
		}
	}
	if (grandmean) {
		list(mu=mu, a=a, d=d, aa=aa, dd=dd)
	} else {
		list(a=a, d=d, aa=aa, dd=dd)
	}
}


fit.2pop.F2 <- function(P1, P2, F1, F2, what=c("a", "d", "aa")) {
	Z <- do.call(cbind, Z.2pop.F2(P1, P2, F1, F2)[c("mu", what)])
	lm(c(P1, P2, F1, F2) ~ 0 + Z)
}


fit.Npop.Fmu <- function(fP1, fP2, fgen, phen, what=c("a", "d", "aa"), grandmean=FALSE) {
	allZ <- Z.Npop.Fmu(fP1, fP2, fgen, grandmean=grandmean)[if (grandmean) c("mu", what) else what]
	Z <- do.call(cbind, allZ)
	colnames(Z) <- do.call(c, lapply(names(allZ), function(x) paste(x, colnames(allZ[[x]]), sep=if(x != "mu") "." else "")))
	if (!grandmean)
		phen <- phen - mean(phen, na.rm=TRUE)
	lm(phen ~ 0 + Z)
}

formula.multi <- function(nn, e.unique=FALSE) {
	# nn is the vector of LINEAR effect names (i.e., including aa and dd)
	mu <- nn[grepl("^mu",    nn)]
	a  <- nn[grepl("^a\\.",  nn)]
	d  <- nn[grepl("^d\\.",  nn)]
	aa <- nn[grepl("^aa\\.", nn)]
	dd <- nn[grepl("^dd\\.", nn)]
	ee.aa <- if (e.unique) "ee" else paste0("ee.", sub("aa.", "", aa))
	ee.dd <- if (e.unique) "ee" else paste0("ee.", sub("dd.", "", dd))
	
	form.mean <- paste0("phen ~ ")
	form.mu   <- paste0("Z[,\"", mu, "\"]*mu")
	form.a    <- paste0( "Z[,\"", a, "\"]*", a, collapse=" + ")
	form.d    <- paste0( "Z[,\"", d, "\"]*", d, collapse=" + ")
	form.aa   <- paste0( "Z[,\"", aa,"\"]*", "a.", sub(".*\\.(.*)x.*", "\\1", aa, perl=TRUE), "*a.", sub(".*x(.*)$", "\\1", aa, perl=TRUE), "*", ee.aa, collapse=" + ")
	form.dd   <- paste0( "Z[,\"", dd,"\"]*", d, "*", d, "*", ee.dd, collapse=" + ")
	
	form <- paste0(
		form.mean, 
		if (length(mu) > 0) form.mu else "0",
		if (length(a)  > 0) paste0(" + ", form.a)  else "", 
		if (length(d)  > 0) paste0(" + ", form.d)  else "", 
		if (length(aa) > 0) paste0(" + ", form.aa) else "", 
		if (length(dd) > 0) paste0(" + ", form.dd) else "")
	form
}

effects.clean <- function(fP1, fP2, fgen, phen, what.lin=c("a","d","aa"), grandmean=FALSE) {
	lin <- fit.Npop.Fmu(fP1, fP2, fgen, phen, what=what.lin, grandmean=grandmean)
	clin <- coef(lin)
	clin <- clin[!is.na(clin)]
	ans <- sub("^Z", "", names(clin))
	ans
}

start.multi <- function(fP1, fP2, fgen, phen, what.multi=c("a", "d", "ee"), e.unique=FALSE, grandmean=FALSE) {
	lin <- fit.Npop.Fmu(fP1, fP2, fgen, phen, what=c("a", "d", "aa")[c("a","d","ee") %in% what.multi], grandmean=grandmean)
	clin <- coef(lin)
	names(clin) <- sub("^Z", "", names(clin))
	
	eff <- effects.clean(fP1, fP2, fgen, phen, what=c("a", "d", "aa")[c("a","d","ee") %in% what.multi], grandmean=grandmean)
	start <- list(
		mu = clin["mu"],
		a  = clin[eff[grep("^a\\.", eff)]],
		d  = clin[eff[grep("^d\\.", eff)]],
		aa = clin[eff[grep("^aa\\.", eff)]])
	if (e.unique) {
		start$ee <- c(ee=0)
	} else {
		start$ee <- setNames(rep(0, length(start$aa)), nm=sub("aa", "ee", names(start$aa)))
	}
	nm <- if(grandmean) c("mu", what.multi) else what.multi
	ans <- setNames(unlist(start[nm]), nm=unlist(sapply(start[nm], names)))
	ans
}

fit.Npop.Fmu.multi <- function(fP1, fP2, fgen, phen, what.lin=c("a", "d", "aa", "dd"), e.unique=FALSE, grandmean=FALSE) {
	what.multi <- c("a", "d", "ee")[c("a" %in% what.lin, "d" %in% what.lin, "aa" %in% what.lin || "dd" %in% what.lin)]
	allZ <- Z.Npop.Fmu(fP1, fP2, fgen, grandmean=grandmean)[if(grandmean) c("mu", what.lin) else what.lin]
	Z <- do.call(cbind, allZ)
	colnames(Z) <- do.call(c, lapply(names(allZ), function(x) paste(x, colnames(allZ[[x]]), sep=if(x == "mu") "" else ".")))
	ss <- start.multi(fP1, fP2, fgen, phen, what.multi, e.unique, grandmean=grandmean)
	form <- formula.multi(effects.clean(fP1, fP2, fgen, phen, what.lin, grandmean=grandmean), e.unique=e.unique)
	if (!grandmean)
		phen <- phen - mean(phen, na.rm=TRUE)
		
	nls(formula=as.formula(form), start=ss, control=nls.control(maxiter=50, warnOnly=TRUE))
}

fit.Npop.Fmu.rand <- function(fP1, fP2, fgen, phen, what=c("a", "d", "aa")) {
	library(hglm)
	# hglm is allergic to NAs
	isna <- is.na(phen)
	fP1  <- fP1 [!isna]
	fP2  <- fP2 [!isna]
	fgen <- fgen[!isna]
	phen <- phen[!isna]
	
	allZ <- Z.Npop.Fmu(fP1, fP2, fgen)[what]   # not mu here
	Z <- do.call(cbind, allZ)
	hglm(X=do.call(cbind, lapply(allZ, rowSums)), y=phen, Z=Z, RandC=sapply(allZ, ncol), calc.like=TRUE)
}


fit.2pop <- function(P1, P2, F1, F2) {
	start <- fit.2pop.noDD(P1, P2, F1, F2)
	ans <- try(optim(par=start, 
	      fn=function(pp) 
	      	(P1 - (pp["zR"] - MODEL["P","a"] * pp["a"])) ^ 2 +
	      	(P2 - (pp["zR"] + MODEL["P","a"] * pp["a"])) ^ 2 +
	      	(F1 - (pp["zR"] + MODEL["F1","d"] * pp["d"] + MODEL["F1","dd"] * pp["e"] * pp["d"]^2  + MODEL["F1","aa"] * pp["e"] * pp["a"]^2 )) ^ 2 +
	      	(F2 - (pp["zR"] + MODEL["F2","d"] * pp["d"] + MODEL["F2","dd"] * pp["e"] * pp["d"]^2  + MODEL["F2","aa"] * pp["e"] * pp["a"]^2 )) ^ 2
	      ))
	if (class(ans) == "try-error") 
		return(rep(NA, 4))
	else
		ans$par
}

fit.2pop.noDD <- function(P1, P2, F1, F2) {
	ans <- c(
		zR = (P2+P1)/2,
		a  = (P2-P1)/2)
	if (ans["a"] == 0) ans["a"] <- 1e-6
	F1s <- (F1-ans["zR"])/ans["a"]
	F2s <- (F2-ans["zR"])/ans["a"]
	ds <- 2*(F1s-F2s)
	es <- 2*(F1s-2*F2s)
	ans["d"] <- ds*ans["a"]
	ans["e"] <- es/ans["a"]
	ans
}
