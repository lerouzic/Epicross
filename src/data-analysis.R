library(parallel)

mc.cores <- min(8, detectCores() - 1)
#~ mc.cores <- 1

source("../src/genet-model.R")
source("../src/admb-helper.R")
source("../src/hglm-helper.R")


dd <- read.table("../data/data_clean.txt", stringsAsFactors=FALSE)

mean.lines <- lapply(setNames(nm=c("Weight", "Fitness")), function(trait) {
	lines <- unique(c(dd$Mother_line, dd$Father_line))
	ans <- list(
		P = sapply(lines, function(line) mean(dd[dd$Mother_line == line & dd$Father_line == line, trait], na.rm=TRUE)), 
		F1 = outer(lines, lines, function(line1, line2) mapply(line1, line2, FUN=function(l1, l2) mean(dd[dd$Mother_line == l1 & dd$Father_line == l2 & dd$Gen == "F1", trait], na.rm=TRUE))),
		F2 = outer(lines, lines, function(line1, line2) mapply(line1, line2, FUN=function(l1, l2) mean(dd[dd$Mother_line == l1 & dd$Father_line == l2 & dd$Gen == "F2", trait], na.rm=TRUE)))
		)
	colnames(ans$F1) <- rownames(ans$F1) <- colnames(ans$F2) <- rownames(ans$F2) <- names(ans$P)
	ans
	})

crosses.lines <- lapply(mean.lines, function(cc) {
	ans <- data.frame()
	for (l1 in names(cc$P))
		for (l2 in names(cc$P)) {
			tmp.ans <- data.frame(P1=cc$P[l1], P2=cc$P[l2], F1=cc$F1[l1,l2], F2=cc$F2[l1,l2], 
								sF1=(cc$F1[l1,l2]-(cc$P[l1]+cc$P[l2])/2)/(abs(cc$P[l2]-cc$P[l1])/2), 
								sF2=(cc$F2[l1,l2]-(cc$P[l1]+cc$P[l2])/2)/(abs(cc$P[l2]-cc$P[l1])/2))
			rownames(tmp.ans) <- paste(l1, l2, sep="x")
			if (l1 != l2 && all(is.finite(unlist(tmp.ans))))
				ans <- rbind(ans, tmp.ans)
		}
	ans
})


mean.pops <- lapply(setNames(nm=c("Weight", "Fitness")), function(trait) {
	pops <- unique(c(dd$Mother_pop, dd$Father_pop))
	ans <- list(
		P = sapply(pops, function(pop) mean(dd[dd$Mother_pop == pop & dd$Father_pop == pop, trait], na.rm=TRUE)), 
		F1 = outer(pops, pops, function(pop1, pop2) mapply(pop1, pop2, FUN=function(l1, l2) mean(dd[dd$Mother_pop == l1 & dd$Father_pop == l2 & dd$Gen == "F1", trait], na.rm=TRUE))),
		F2 = outer(pops, pops, function(pop1, pop2) mapply(pop1, pop2, FUN=function(l1, l2) mean(dd[dd$Mother_pop == l1 & dd$Father_pop == l2 & dd$Gen == "F2", trait], na.rm=TRUE)))
		)
	colnames(ans$F1) <- rownames(ans$F1) <- colnames(ans$F2) <- rownames(ans$F2) <- names(ans$P)
	ans
	})

crosses.pops <- lapply(mean.pops, function(cc) {
	ans <- data.frame()
	for (l1 in names(cc$P))
		for (l2 in names(cc$P)) {
			tmp.ans <- data.frame(P1=cc$P[l1], P2=cc$P[l2], F1=cc$F1[l1,l2], F2=cc$F2[l1,l2], 
								sF1=(cc$F1[l1,l2]-(cc$P[l1]+cc$P[l2])/2)/(abs(cc$P[l2]-cc$P[l1])/2), 
								sF2=(cc$F2[l1,l2]-(cc$P[l1]+cc$P[l2])/2)/(abs(cc$P[l2]-cc$P[l1])/2))
			rownames(tmp.ans) <- paste(l1, l2, sep="x")
			if (l1 != l2 && all(is.finite(unlist(tmp.ans))))
				ans <- rbind(ans, tmp.ans)
		}
	ans
})

col <- c(
		ref="darkgray", 
		F1 ="orange", 
		F2 ="green",
		D  ="seagreen3",
		AxA="tomato",
		add="orchid2")

crossfig1 <- function(dat, lwd=3, pch=1, ...) {
	plot(0.5*(dat$P1 + dat$P2), dat$F1, pch=pch, col=col["F1"], xlab="Mid-parent", ylab="Offspring", ylim=range(c(dat$F1, dat$F2)), ...)
	points(0.5*(dat$P1 + dat$P2), dat$F2, pch=pch, col=col["F2"])
	abline(a=0, b=1, col=col["ref"])
	abline(lm(dat$F1 ~ I(0.5*(dat$P1 + dat$P2))), col=col["F1"], lwd=lwd)
	abline(lm(dat$F2 ~ I(0.5*(dat$P1 + dat$P2))), col=col["F2"], lwd=lwd)
	abline(v=mean(0.5*(dat$P1 + dat$P2)), lty=2, col=col["ref"])
	legend("topright", lty=1, col=col[c("F1","F2","ref")], lwd=c(lwd, lwd, 1), legend=c("F1","F2","Additive"))
}

crossfig2 <- function(dat, lwd=3, pch=1, ...) {
	plot(dat$F1 - 0.5*(dat$P1 + dat$P2), dat$F1 - dat$F2, pch=pch,  xlab="F1 - Parents = 2D-AA", ylab="F1 - F2 = D",  ...)
	abline(a=0, b=0.5, col=col["D"], lwd=2)
	abline(a=0, b=0, lwd=2, col=col["AxA"])
	points(0, 0, pch=16, col=col["add"], cex=3)
	points(mean(dat$F1 - 0.5*(dat$P1 + dat$P2)), mean(dat$F1 - dat$F2), cex=3, pch=1, lwd=3)
	legend("topleft", lty=c(0,0,1,1), pch=c(1, 16, NA, NA), col=c("black", col["add"], col["D"], col["AxA"]), legend=c("Data", "Additive", "Dominance", "AxA epistasis"))
}

crossfig3 <- function(dat, lwd=3, pch=1, ...) {
	plot(dat$F1 -2*dat$F2 + 0.5*(dat$P1 + dat$P2), dat$F1 - dat$F2, pch=pch,  xlab="Parents + F1 - 2 F2 = AA", ylab="F1 - F2 = D",  ...)
	abline(v=0, col=col["D"], lwd=2)
	abline(h=0, col=col["AxA"], lwd=2)
	points(0, 0, pch=16, col=col["add"], cex=3)
	points(mean(dat$F1 -2*dat$F2 + 0.5*(dat$P1 + dat$P2)), mean(dat$F1 - dat$F2), cex=3, pch=1, lwd=3)
	legend("topleft", lty=c(0,0,1,1), pch=c(1, 16, NA, NA), col=c("black", col["add"], col["D"], col["AxA"]), legend=c("Data", "Additive", "Dominance", "AxA epistasis"))
}


crossfig2.theor <- function( ...) {
	plot(NULL, xlim=c(-1, 1), ylim=c(-1, 1), xaxt="n", yaxt="n", xlab="F1 - Parents = 2D-AA", ylab="F1 - F2 = D",  ...)
	axis(1, at=0)
	axis(2, at=0)
	abline(a=0, b=0.5, col=col["D"], lwd=2)
	abline(a=0, b=0, lwd=2, col=col["AxA"])
	points(0, 0, pch=16, col=col["add"], cex=3)
	
	text(-0.5, 0, "No dominance", pos=3, col=col["AxA"])
	text(0.5, 0.25, "No epistasis", pos=3, col=col["D"], srt=45/2)
	text(0.25, -0.25, "Additive", col=col["add"])
	arrows(x0=0.15, x1=0.05, y0=-0.2, y1=-0.05, col=col["add"], length=0.05)
	
	text(0.2, -0.5, "D < 0\nAxA < 0", col=col["ref"])
	text(-0.2, 0.5, "D > 0\nAxA > 0", col=col["ref"])
	text(0.7, 0.15, "D > 0\nAxA < 0", col=col["ref"])
	text(-0.7, -0.15, "D < 0\nAxA > 0", col=col["ref"])
}

pdf("../results/DataF1lines.pdf", width=12, height=6)
	layout(t(1:2))
	crossfig1(crosses.lines$Weight, main="Weight (lines)")
	crossfig1(crosses.lines$Fitness, main="Siliques (lines)")
dev.off()

pdf("../results/DataF1pops.pdf", width=12, height=6)
	layout(t(1:2))
	crossfig1(crosses.pops$Weight, main="Weight (populations)")
	crossfig1(crosses.pops$Fitness, main="Siliques (populations)")
dev.off()

pdf("../results/DataF2lines.pdf", width=12, height=6)
	layout(t(1:2))
	crossfig3(crosses.lines$Weight, main="Weight (lines)")
	crossfig3(crosses.lines$Fitness, main="Siliques (lines)")
dev.off()

pdf("../results/DataF2pops.pdf", width=12, height=6)
	layout(t(1:2))
	crossfig3(crosses.pops$Weight, main="Weight (populations)")
	crossfig3(crosses.pops$Fitness, main="Siliques (populations)")
dev.off()

pdf("../results/DataF2theor.pdf", width=6, height=6)
	crossfig2.theor()
dev.off()
