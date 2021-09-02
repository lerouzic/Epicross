
library(hglm)
source("../src/hglm-copy.R") # Fixing minor issues in hglm

source("../src/genet-model.R")

.check.par <- function(pop.par1, pop.par2, gen, effect) {
	stopifnot(
		gen %in% c("F1","F2"),
		is.vector(pop.par1), is.vector(pop.par2),
		length(pop.par1) == length(pop.par2),
		effect %in% colnames(MODEL))
}

structure.matrix.1D <- function(pop.par1, pop.par2, gen, effect="a") {
	.check.par(pop.par1, pop.par2, gen, effect)
	
	lv <- gtools::mixedsort(levels(factor(c(as.character(pop.par1), as.character(pop.par2)))))
	
	X <- matrix(0, nrow=length(pop.par1), ncol=length(lv))
	colnames(X) <- lv
	gen <- ifelse(pop.par1 == pop.par2, "P", gen)
	
	for (i in seq_along(pop.par1))
		X[i, as.character(pop.par1[i])] <- X[i, as.character(pop.par1[i])] + MODEL[gen[i], effect]
	for (i in seq_along(pop.par2))
		X[i, as.character(pop.par2[i])] <- X[i, as.character(pop.par2[i])] + MODEL[gen[i], effect]
	# Removing empty columns
	X[,colSums(abs(X))>0] 
}

structure.matrix.2D <- function(pop.par1, pop.par2, gen, effect="d") {
	.check.par(pop.par1, pop.par2, gen, effect)
	
	# Not very efficient algorithm. 
	lv <- gtools::mixedsort(levels(factor(c(as.character(pop.par1), as.character(pop.par2)))))
	lv.inter <- outer(lv, lv, paste, sep="-")
	lv.inter <- c(lv.inter[upper.tri(lv.inter, diag=TRUE)])
	
	X <- matrix(0, nrow=length(pop.par1), ncol=length(lv.inter))
	colnames(X) <- lv.inter
	gen <- ifelse(pop.par1 == pop.par2, "P", gen)
	
	for (i in seq_along(pop.par1)) {
		nm <- paste(as.character(pop.par1[i]), as.character(pop.par2[i]), sep="-")
		if (!nm %in% colnames(X))
			nm <- paste(as.character(pop.par2[i]), as.character(pop.par1[i]), sep="-")
		X[i, nm] <- X[i, nm] + MODEL[gen[i], effect]
	}
	# Removing empty columns
	X[,colSums(abs(X))>0]
}

additive.effects <- function(pop.par1, pop.par2, gen) {
	list(
		X=NULL,
		Z=structure.matrix.1D(pop.par1, pop.par2, gen, "a")
	)
}

dominance.effects <- function(pop.par1, pop.par2, gen) {
	Zdom <- structure.matrix.2D(pop.par1, pop.par2, gen, "d")
	list(
		X=as.matrix(rowSums(Zdom)),
		Z=Zdom
	)
}

AxA.effects <- function(pop.par1, pop.par2, gen) {
	Zaa <- structure.matrix.2D(pop.par1, pop.par2, gen, "aa")
	list(
		X=as.matrix(rowSums(Zaa)),
		Z=Zaa
	)
}

DxD.effects <- function(pop.par1, pop.par2, gen) {
	Zdd <- structure.matrix.2D(pop.par1, pop.par2, gen, "dd")
	list(
		X=as.matrix(rowSums(Zdd)),
		Z=Zdd
	)
}

hglm.wrapper <- function(dd, trait="Weight", A.pop=TRUE, A.line=FALSE, D.pop=TRUE, D.line=FALSE, E.pop=FALSE, E.line=FALSE, AA.pop=TRUE, AA.line=FALSE, DD.pop=TRUE, DD.line=FALSE) {
	
	if (!E.pop)  { AA.pop  <- DD.pop  <- FALSE }
	if (!E.line) { AA.line <- DD.line <- FALSE }
	if (E.pop  && A.pop)  AA.pop <- TRUE
	if (E.pop  && D.pop)  DD.pop <- TRUE
	if (E.line && A.line) AA.line <- TRUE
	if (E.line && D.line) DD.line <- TRUE
	
	X <- matrix(rep(1, nrow(dd)))
	
	listX <- list(mu=X)
	listZ <- list()
		
	if (A.pop) {
		ssA <- additive.effects(dd$Mother_pop, dd$Father_pop, dd$Gen)
		listZ <- c(listZ, list(pop.a=ssA$Z))
	}
	
	if (A.line) {
		ssA <- additive.effects(dd$Mother_line, dd$Father_line, dd$Gen)
		listZ <- c(listZ, list(lin.a=ssA$Z))
	}
	
	if (D.pop) {
		ssD <- dominance.effects(dd$Mother_pop, dd$Father_pop, dd$Gen)
		listX <- c(listX, list(pop.d=ssD$X))
		listZ <- c(listZ, list(pop.d=ssD$Z))
	}
	
	if (D.line) {
		ssD <- dominance.effects(dd$Mother_line, dd$Father_line, dd$Gen)
		listX <- c(listX, list(lin.d=ssD$X))
		listZ <- c(listZ, list(lin.d=ssD$Z))
	}
	
	if (AA.pop) {
		ssAA <- AxA.effects(dd$Mother_pop, dd$Father_pop, dd$Gen)
		listX <- c(listX, list(pop.aa=ssAA$X))
		listZ <- c(listZ, list(pop.aa=ssAA$Z))
	}
	
	if (AA.line) {
		ssAA <- AxA.effects(dd$Mother_line, dd$Father_line, dd$Gen)
		listX <- c(listX, list(lin.aa=ssAA$X))
		listZ <- c(listZ, list(lin.aa=ssAA$Z))
	}
	
	if (DD.pop) {
		ssDD <- DxD.effects(dd$Mother_pop, dd$Father_pop, dd$Gen)
		listX <- c(listX, list(pop.dd=ssDD$X))
		listZ <- c(listZ, list(pop.dd=ssDD$Z))
	}
	
	if (DD.line) {
		ssDD <- DxD.effects(dd$Mother_line, dd$Father_line, dd$Gen)
		listX <- c(listX, list(lin.dd=ssDD$X))
		listZ <- c(listZ, list(lin.dd=ssDD$Z))
	}
	
	y <- dd[,trait]
	NAy <- is.na(y)
	y <- y[!NAy]
	
	X <- do.call(cbind, listX)
	colnames(X) <- unlist(lapply(names(listX), function(nn) { cc <- colnames(listX[[nn]]); if (is.null(cc)) nn else paste(nn, cc, sep=".")}))
	X <- X[!NAy,,drop=FALSE]
	
	Z <- do.call(cbind, listZ)
	colnames(Z) <- unlist(lapply(names(listZ), function(nn) { cc <- colnames(listZ[[nn]]); if (is.null(cc)) nn else paste(nn, cc, sep=".")}))
	Z <- Z[!NAy,]
	
	startval <- c(
		c(mean(y), rep(0, ncol(X)-1)),
		rep(0, ncol(Z)),
		rep(1, length(listZ)),
		var(y))

	hglm(y=y, X=X, Z=Z, RandC=sapply(listZ, ncol), startval=startval, verbose=FALSE)
}

hglm.summ <- function(obj, what="estimates") {
	stopifnot("hglm" %in% class(obj))
	
	if (what == "estimates") 
		ans <- list(
			mu     = obj$fixef["mu"],
			pop.d  = obj$fixef["pop.d"],
			pop.aa = obj$fixef["pop.aa"],
			pop.dd = obj$fixef["pop.dd"],
			lin.d  = obj$fixef["lin.d"],
			lin.aa = obj$fixef["lin.aa"],
			lin.dd = obj$fixef["lin.dd"],
			Vpop.a = obj$varRanef[which(names(obj$RandC) == "pop.a")],
			Vlin.a = obj$varRanef[which(names(obj$RandC) == "lin.a")],
			Vpop.d = obj$varRanef[which(names(obj$RandC) == "pop.d")],
			Vlin.d = obj$varRanef[which(names(obj$RandC) == "lin.d")],
			Vpop.aa= obj$varRanef[which(names(obj$RandC) == "pop.aa")],
			Vlin.aa= obj$varRanef[which(names(obj$RandC) == "lin.aa")],
			Vpop.dd= obj$varRanef[which(names(obj$RandC) == "pop.dd")],
			Vlin.dd= obj$varRanef[which(names(obj$RandC) == "lin.dd")],
			VR     = obj$varFix
		)
	if (what == "variances")
		ans <- list(
		)
	lapply(ans, function(x) if(length(x) == 0) NA else x)
}

