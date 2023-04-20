source("../src/common.R")


dd <- read.table(data.file, stringsAsFactors=FALSE)

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


crossfig2 <- function(dat, lwd=3, pch=1, ...) {
	plot(dat$F1 - 0.5*(dat$P1 + dat$P2), dat$F1 - dat$F2, pch=pch,  xlab="F1 - Parents = 2D-AA", ylab="F1 - F2 = D",  ...)
	abline(a=0, b=0.5, col=col["D"], lwd=2)
	abline(a=0, b=0, lwd=2, col=col["AxA"])
	points(0, 0, pch=16, col=col["add"], cex=3)
	points(mean(dat$F1 - 0.5*(dat$P1 + dat$P2)), mean(dat$F1 - dat$F2), cex=3, pch=1, lwd=3)
	legend("topleft", lty=c(0,0,1,1), pch=c(1, 16, NA, NA), col=c("black", col["add"], col["D"], col["AxA"]), legend=c("Data", "Additive", "Dominance", "AxA epistasis"))
}




pdf("../results/Fig1BC.pdf", width=12, height=6)
	layout(t(1:2))
	crossfig2(crosses.lines$Weight, main="Weight")
	crossfig2(crosses.lines$Fitness, main="Siliques")
dev.off()


pdf("../results/FigS1.pdf", width=12, height=6)
	layout(t(1:2))
	crossfig2(crosses.pops$Weight, main="Weight")
	crossfig2(crosses.pops$Fitness, main="Siliques")
dev.off()
