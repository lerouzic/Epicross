
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

Z.Npop.Fmu <- function(fP1, fP2, fgen) {
	stopifnot(length(fP1) == length(fP2), length(fP1) == length(fgen))
	
	n <- length(fP1)                # Number of measurements
	f <- factor(unique(c(fP1,fP2))) # Populations
	i <- factor(unique(paste(as.character(fP1), as.character(fP2), sep="-")[fP1 != fP2]))
	
	mu <- rep(1, length(fP1))
	a  <- matrix(0, nrow=n, ncol=length(f))
	d  <- matrix(0, nrow=n, ncol=length(i))
	aa <- matrix(0, nrow=n, ncol=length(i))
	dd <- matrix(0, nrow=n, ncol=length(i))
	
	for (j in seq_along(fP1)) {
		fi1 <- which(f == fP1[j])
		fi2 <- which(f == fP2[j])
		
		a[j, fi1] <- a[j, fi1] + 1/2
		a[j, fi2] <- a[j, fi2] + 1/2
		
		if (fP1[j] != fP2[j]) {
			ii <- which(paste(fP1[j], fP2[j], sep="-") %in% i)
			aa[j,ii] <- -1
			if (fgen[j] == "F1") {
				d[j,ii]  <- 2
			} else {
				d[j,ii]  <- 1
				dd[j,ii] <- -1
			}
		}
	}
	list(mu=mu, a=a, d=d, aa=aa, dd=dd)
}


fit.2pop.F2 <- function(P1, P2, F1, F2, what=c("a", "d", "aa")) {
	Z <- do.call(cbind, Z.2pop.F2(P1, P2, F1, F2)[c("mu", what)])
	lm(c(P1, P2, F1, F2) ~ 0 + Z)
}


fit.Npop.Fmu <- function(fP1, fP2, fgen, phen, what=c("a", "d", "aa") ) {
	Z <- do.call(cbind, Z.Npop.Fmu(fP1, fP2, fgen)[c("mu", what)])
	lm(phen ~ 0 + Z)
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
