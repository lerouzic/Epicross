source("../src/gp-simulator.R")
source("../src/genet-model.R")


scalecross <- function(datacross) {
	toswap <- datacross$P1 > datacross$P2
	datacross$P1.id <- ifelse(toswap, datacross$P2.id, datacross$P1.id)
	datacross$P2.id <- ifelse(toswap, datacross$P1.id, datacross$P2.id)
	datacross$P1 <- ifelse(toswap, datacross$P2, datacross$P1)
	datacross$P2 <- ifelse(toswap, datacross$P1, datacross$P2)
	
	datacross$sF1 <- (datacross$F1 - (datacross$P1+datacross$P2)/2)/((datacross$P1-datacross$P2)/2)
	datacross$sF2 <- (datacross$F2 - (datacross$P1+datacross$P2)/2)/((datacross$P1-datacross$P2)/2)
	datacross$sF1[!is.finite(datacross$sF1)] <- NA
	datacross$sF2[!is.finite(datacross$sF2)] <- NA
	datacross
}

plotscross <- function(datacross, xlim=NULL, ylim=NULL, ...) {
	datacross <- scalecross(datacross)
	if (is.null(xlim)) xlim <- range(datacross$sF1, na.rm=TRUE)
	if (is.null(ylim)) ylim <- range(datacross$sF2, na.rm=TRUE)
	plot(NULL, xlim=xlim, ylim=ylim, xlab="F1 (relative to mid-parent)", ylab="F2 (relative to mid-parent)", asp=1, ...)
	points(datacross$sF1, datacross$sF2)
	abline(v=0, h=0, col="gray")
	abline(a=0, b=1, col="darkgreen", lwd=2)
	if (mean(datacross$sF1, na.rm=TRUE) > 0)
		capt.loc <- 0.7*max(xlim)
	else
		capt.loc <- 0.7*min(xlim)
	text(capt.loc,capt.loc,"No dominance", pos=3, srt=45, col="darkgreen")
	abline(a=0, b=0.5, col="darkviolet", lwd=2)
	text(capt.loc,capt.loc/2,"No epistasis", pos=3, srt=28, col="darkviolet")
}

plotcross2 <- function(datacross, xlim=NULL, ylim=NULL, jitter=0.005, maxepi=10, epsilon=FALSE, ...) {
	effects <- t(mapply(datacross$P1, datacross$P2, datacross$F1, datacross$F2, FUN=fit.2pop))
	
	if (epsilon) { 
		effects[effects[,"e"] < -maxepi ,"e"] <- -maxepi
		effects[effects[,"e"] >  maxepi ,"e"] <-  maxepi
	} else {
		effects[,"e"] <- effects[,"e"] * effects[,"a"]^2
	}
	if (is.null(xlim)) xlim <- c(min(min(effects[,"d"]), -0.1), max(max(effects[,"d"]), 0.1))
	if (is.null(ylim)) ylim <- c(min(min(effects[,"e"]), -0.1), max(max(effects[,"e"]), 0.1))

	plot(effects[,"d"]+rnorm(nrow(effects), 0, sd=jitter*diff(xlim)), effects[,"e"]+rnorm(nrow(effects), 0, sd=jitter*diff(ylim)), xlim=xlim, ylim=ylim, xlab="Dominance D", ylab=if (epsilon) expression("Epistasis "*epsilon) else expression ("Scaled epistasis a"^2*epsilon), ...)
	abline(h=0, v=0, col="darkgray")
	if (epsilon) {
		if (any(effects[,"e"] == maxepi))
			abline(h=maxepi, col="darkgray", lty=2)
		if (any(effects[,"e"] == -maxepi))
			abline(h=-maxepi, col="darkgray", lty=2)
	}
}



nloc <- 6
npop <- 12

simA <- meancross.multilin(npop, nloc, a=rep(1, nloc))
plotcross2(simA, main = paste0("Additive, ", nloc, " loci"))


simD0 <- meancross.multilin(npop, nloc, a=rep(1, nloc), d=seq(-1, 1, length=nloc))
plotcross2(simD0, main = paste0("Zero-average dominance, ", nloc, " loci"))

simDP <- meancross.multilin(npop, nloc, a=rep(1, nloc), d=rep(1, nloc))
plotcross2(simDP, main = paste0("Positive constant dominance, ", nloc, " loci"))

simDP2 <- meancross.multilin(npop, nloc, a=rep(1, nloc), d=seq(0, 2, length=nloc))
plotcross2(simDP2, main = paste0("Positive variable dominance, ", nloc, " loci"))


simE0 <- meancross.multilin(npop, nloc, a=rep(1, nloc), d=rep(0, nloc), e=rnorm(nloc*(nloc-1)/2, 0, 0.1))
plotcross2(simE0, main = paste0("Zero-average epistasis, ", nloc, " loci"))

simEP <- meancross.multilin(npop, nloc, a=rep(1, nloc), d=rep(0, nloc), e=rep(1, nloc*(nloc-1)/2))
plotcross2(simEP, main = paste0("Positive constant epistasis, ", nloc, " loci"))
