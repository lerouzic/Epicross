
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

fit.Npop.Fmu <- function(fP1, fP2, fgen, phen, what=c("a", "d", "aa"), grandmean=FALSE) {
	allZ <- Z.Npop.Fmu(fP1, fP2, fgen, grandmean=grandmean)[if (grandmean) c("mu", what) else what]
	Z <- do.call(cbind, allZ)
	colnames(Z) <- do.call(c, lapply(names(allZ), function(x) paste(x, colnames(allZ[[x]]), sep=if(x != "mu") "." else "")))
	if (!grandmean)
		phen <- phen - mean(phen, na.rm=TRUE)
	lm(phen ~ 0 + Z)
}

filter.effects <- function(model, what="a", remove.Zname=FALSE) {
	mm <- coef(model)
	ans <- mm[grep(paste0("^Z?", what, "\\."), names(mm))]
	if (remove.Zname)
		names(ans) <- sapply(strsplit(names(ans), split="\\."), "[", 2)
	ans
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
