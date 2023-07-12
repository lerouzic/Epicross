source("../src/common.R")

sep.lines <- 0.4

crossfig2.theor <- function( ...) {
	plot(NULL, xlim=c(-1, 1), ylim=c(-1, 1), xaxt="n", yaxt="n", xlab="F1 - Parents = 2D-AA", ylab="F1 - F2 = D",  ...)
	axis(1, at=0)
	axis(2, at=0)
	abline(a=0, b=0.5, col=col["D"], lwd=2)
	abline(a=0, b=0, lwd=2, col=col["AxA"])
	points(0, 0, pch=16, col=col["add"], cex=3)
	
	text(-0.5, 0, "No dominance (D=0)", pos=3, col=col["AxA"])
	text(0.5, 0.25, "No epistasis (AA = 0)", pos=3, col=col["D"], srt=45/2)
	text(0.25, -0.25, "Additive", col=col["add"])
	arrows(x0=0.15, x1=0.05, y0=-0.2, y1=-0.05, col=col["add"], length=0.05)
	
	text(0.2, -0.5, "D < 0\nAA < 0", col=col["ref"])
	text(-0.2, 0.5, "D > 0\nAA > 0", col=col["ref"])
	text(0.7, 0.15, "D > 0\nAA < 0", col=col["ref"])
	text(-0.7, -0.15, "D < 0\nAA > 0", col=col["ref"])
	
	for (d in seq(sep.lines,2,by=sep.lines)) {
		abline(h=c(d, -d), lty=3, col=adjustcolor(col["D"], 0.9))
	}
	
	for (aa in seq(sep.lines,2,by=sep.lines)) {
		abline(a=aa, b=0.5, lty=3, col=adjustcolor(col["AxA"], 0.9))
		abline(a=-aa, b=0.5, lty=3, col=adjustcolor(col["AxA"], 0.9))
	}
	
}

pdf("../results/Fig1A.pdf", width=6, height=6)
	crossfig2.theor()
dev.off()
