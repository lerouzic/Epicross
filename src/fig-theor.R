source("../src/common.R")

sep.lines <- 0.4

crossfig2.theor <- function(cex.text=0.8, ...) {
	plot(NULL, xlim=c(-1, 1), ylim=c(-1, 1), xaxt="n", yaxt="n", xlab="F1 - Parents = 2D-AA", ylab="F1 - F2 = D",  ...)
	axis(1, at=0)
	axis(2, at=0)
	abline(a=0, b=0.5, col=col["D"], lwd=2)
	abline(a=0, b=0, lwd=2, col=col["AxA"])
	points(0, 0, pch=16, col=col["add"], cex=3)
	
	text(-0.5, 0, "No dominance (D=0)", pos=3, col=col["AxA"], cex=cex.text)
	text(0.5, 0.25, "No epistasis (AA = 0)", pos=3, col=col["D"], srt=45/2, cex=cex.text)
	text(0.25, -0.25, "Additive", col=col["add"], cex=cex.text)
	arrows(x0=0.15, x1=0.05, y0=-0.2, y1=-0.05, col=col["add"], length=0.05)
	
	text(0.2, -0.5, "D < 0\nAA < 0", col=col["ref"], cex=cex.text)
	text(-0.2, 0.5, "D > 0\nAA > 0", col=col["ref"], cex=cex.text)
	text(0.7, 0.15, "D > 0\nAA < 0", col=col["ref"], cex=cex.text)
	text(-0.7, -0.15, "D < 0\nAA > 0", col=col["ref"], cex=cex.text)
	
	for (d in seq(sep.lines,2,by=sep.lines)) {
		abline(h=c(d, -d), lty=3, col=adjustcolor(col["D"], 0.9))
	}
	
	for (aa in seq(sep.lines,2,by=sep.lines)) {
		abline(a=aa, b=0.5, lty=3, col=adjustcolor(col["AxA"], 0.9))
		abline(a=-aa, b=0.5, lty=3, col=adjustcolor(col["AxA"], 0.9))
	}
	
}

pdf("../results/Fig1A.pdf", width=fig.width/3, height=fig.height, pointsize=fig.pointsize)
	par(mar=fig.mar)
	crossfig2.theor()
dev.off()

simdyn.2loc <- function(mu=0, a=1, d=0, aa=0, self=0, g=50, sel=0.2, start="F2") {
	# sel is the fitness difference between the largest and the lowest phenotype
	
	genot.freq <- function(p1, p2=p1) {
		# Equilibrium genotype frequencies with selfing
		f1 <- c(p1*(1-p1)*self/(2-self) + p1^2, 2*p1*(1-p1)*(1-self)/(1-self/2), p1*(1-p1)*self/(2-self)+(1-p1)^2)
		f2 <- c(p2*(1-p2)*self/(2-self) + p2^2, 2*p2*(1-p2)*(1-self)/(1-self/2), p2*(1-p2)*self/(2-self)+(1-p2)^2)
		outer(f1, f2)
	}
	
	# Genotype-Phenotype map
	G <- mu + rbind(
		c(-a, -a/2 + d - aa, - 2*aa), 
		c(-a/2 + d - aa, 2*d - aa, a/2 + d - aa), 
		c(-2*aa, a/2 + d - aa, a))
	# This is based on the traditional 1-locus F2 model
	# (aa: -a - d/2; Aa: +d/2, AA: a - d/2),
	# When expanded to 2 loci, aditive effects were divided by 2 
	# This map is then shifted by  + d - aa (which does not change the meaning of the effects),
	# so that parental populations (double homozygotes) are -a and +a. The reference mu is thus the mid-parent,
	# but coefficients keep the same meaning as in the F2 model.
		
	# Fitness
	W <- (G - min(G))/diff(range(G)) # Between 0 and 1
	W <- 1 + (W-1)*sel
	
	# Genotype frequencies
	F <- if(start == "F2") genot.freq(1/2, 1/2)
	
	ans <- c(sum(F*G), rep(NA, g))
	
	for (t in 1:g) {
		F <- F*W / sum(F*W)                                     # selection
		F <- genot.freq(sum(colSums(F)*c(1,1/2,0)), sum(rowSums(F)*c(1,1/2,0))) # reproduction
		ans[t+1] <- sum(F*G)
	}
	ans
}

pdf("../results/Fig3.pdf", width=fig.width, height=fig.height, pointsize=fig.pointsize)
	gg <- 50 # number of generations
	
	layout(t(c(1,2)))
	par(mar=fig.mar, cex=1)

	self <- c(black=0.95, darkgray=0)
	
	plot(NULL, xlim=c(0,gg), ylim=c(0.42, 0.55), xlab="Generations", ylab=weight.name)
	lines(simdyn.2loc(mu=0.43, a=0.11, d=0, aa=0, self=self[1], g=gg), col=names(self)[1], lty=2)
	for (ns in rev(names(self)))
		lines(simdyn.2loc(mu=0.43, a=0.11, d=0.03, aa=0, self=self[ns], g=gg), col=ns)
	legend("bottomright", lty=c(rep(1,length(self)), 2), col=c(names(self), names(self)[1]), legend=c(paste0("A + D, selfing=", self), paste0("A only, selfing=", self[1])))


	plot(NULL, xlim=c(0,gg), ylim=c(650, 1100), xlab="Generations", ylab=silique.name)
	lines(simdyn.2loc(mu=779, a=262, d=0, aa=0, self=self[1], g=gg), col=names(self)[1], lty=2)
	for (ns in rev(names(self)))
		lines(simdyn.2loc(mu=779, a=262, d=0, aa=83, self=self[ns], g=gg), col=ns)
	legend("bottomright", lty=c(rep(1,length(self)), 2), col=c(names(self), names(self)[1]), legend=c(paste0("A + AA, selfing=", self), paste0("A only, selfing=", self[1])))

dev.off()
