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

sim1 <- meancross.multilin(8, 6, a=1:6)
plotscross(sim1)
mapply(sim1$P1, sim1$P2, sim1$F1, sim1$F2, FUN=fit.2pop)

sim2 <- meancross.multilin(8, 6, a=1:6, d=rnorm(6, 0, 1))
plotscross(sim2, main="Dominance only")
mapply(sim2$P1, sim2$P2, sim2$F1, sim2$F2, FUN=fit.2pop)


sim2b <- meancross.multilin(8, 6, a=rnorm(6, 10, 0.1), d=rnorm(6, 5, 0.1))
plotscross(sim2b)
md <- mapply(sim2b$P1, sim2b$P2, sim2b$F1, sim2b$F2, FUN=fit.2pop)


sim3 <- meancross.multilin(8, 6, a=1:6, epsilon=rnorm(15, 0, 1))
plotscross(sim3, main="Epistasis only")

sim3b <- meancross.multilin(8, 6, a=rnorm(6,0,0.1), epsilon=rnorm(15, 1, 0.))
plotscross(sim3b, xlim=c(-1,1), ylim=c(-1,1))


sim4 <- meancross.multilin(8, 6, a=1:6, d=rnorm(6, 0, 1), epsilon=rnorm(15, 0, 1))
plotscross(sim4)
