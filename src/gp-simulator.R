# Functions to simulate epistatic GP maps

library(noia)

gpmap.multilin <- function(nloc=2, a=rnorm(nloc, 1/nloc, 0), d=rnorm(nloc, 0, 0), epsilon=rnorm(nloc*(nloc-1)/2, 0, 0), AxA=TRUE, AxD=TRUE, DxD=TRUE, reference="Finf") {
	E <- setNames(rep(0, 3^nloc), nm=noia:::effectsNamesGeneral(nloc))
	E[sapply(1:nloc, function(loc) noia:::effNames("a",loc,nloc))] <- a
	E[sapply(1:nloc, function(loc) noia:::effNames("d",loc,nloc))] <- d
	combs <- combn(1:nloc, 2)
	if (AxA)
		E[sapply(1:ncol(combs), function(cc) noia:::effNames(c("a","a"), combs[,cc], nloc))] <- a[combs[1,]] * a[combs[2,]] * epsilon
	if (AxD) {
		E[sapply(1:ncol(combs), function(cc) noia:::effNames(c("a","d"), combs[,cc], nloc))] <- a[combs[1,]] * d[combs[2,]] * epsilon
		E[sapply(1:ncol(combs), function(cc) noia:::effNames(c("d","a"), combs[,cc], nloc))] <- d[combs[1,]] * a[combs[2,]] * epsilon
	}
	if (DxD)
		E[sapply(1:ncol(combs), function(cc) noia:::effNames(c("d","d"), combs[,cc], nloc))] <- d[combs[1,]] * d[combs[2,]] * epsilon
	Smat <- Reduce("%x%", replicate(nloc, noia:::Sloc(reference), simplify=FALSE), init=1)
	if (reference == "Finf") # Finf in noia was not properly implemented for multilocus GP maps
		Smat <- t(t(Smat) - c(0, (0.5*Smat[1,] + 0.5*Smat[nrow(Smat),])[-1]))
	GP <- c(Smat %*% E)
	names(GP) <- as.character(noia:::genNames(nloc))
	GP
}

cross.gen  <- function(genot1, genot2) {
	g1 <- strsplit(genot1, "")
	g2 <- strsplit(genot2, "")
	stopifnot(all(sapply(g1, length) == sapply(g2, length)), all(c(unlist(g1), unlist(g2)) != "2"))
	mapply(g1, g2, FUN=function(gg1, gg2) paste(ifelse(gg1 == gg2, gg1, "2"), collapse=""))
}

cross.freq <- function(genot1, genot2) {
	prob1 <- matrix(c(1,1/2,0,1/2,1/4,0  ,0,0  ,0), ncol=3, dimnames=list(c("1","2","3"), c("1","2","3")))
	prob2 <- matrix(c(0,1/2,1,1/2,1/2,1/2,1,1/2,0), ncol=3, dimnames=list(c("1","2","3"), c("1","2","3")))
	prob3 <- matrix(c(0,0  ,0,0  ,1/4,1/2,0,1/2,1), ncol=3, dimnames=list(c("1","2","3"), c("1","2","3")))
	
	g1 <- unlist(strsplit(genot1, ""))
	g2 <- unlist(strsplit(genot2, ""))
	stopifnot(length(g1) == length(g2))
	gn <- do.call(rbind, strsplit(as.character(noia:::genNames(length(g1))), split=""))
	gcross <- Reduce("%x%", rev(lapply(1:length(g1), function(i) c(prob1[g1[i],g2[i]], prob2[g1[i],g2[i]], prob3[g1[i],g2[i]]))), init=1)
	names(gcross) <- as.character(noia:::genNames(length(g1)))
	gcross
}

meancross.multilin <- function(nlines=10, nloc=4, force.nloc=FALSE, ...) { # ... are arguments for gmap.multilin
	if (nloc > 8 && !force.nloc) 
		stop("Too many loci for a comfortable run, use force.nloc=TRUE to overcome.")
	stopifnot(nlines <= 2^nloc)
	
	gg <- gpmap.multilin(nloc=nloc, ...)

	genot.lines <- sample(names(gg)[!grepl("2", names(gg))], nlines)
	
	comb <- matrix(genot.lines[combn(1:nlines, 2)], nrow=2)
	data.frame(
		P1.id = comb[1,],
		P2.id = comb[2,],
		P1 = gg[comb[1,]], 
		P2 = gg[comb[2,]],
		F1 = gg[cross.gen(comb[1,], comb[2,])],
		F2 = apply(comb, 2, function(ll) sum(gg * cross.freq(cross.gen(ll[1],ll[2]),cross.gen(ll[1],ll[2])))))
}
