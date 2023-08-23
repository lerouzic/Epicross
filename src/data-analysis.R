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
	x <- dat$F1 - 0.5*(dat$P1 + dat$P2)
	y <- dat$F1 - dat$F2
	plot(x, y, pch=pch,  xlab="F1 - Parents = 2D-AA", ylab="F1 - F2 = D",  ...)
	abline(a=0, b=0.5, col=col["D"], lwd=2)
	abline(a=0, b=0, lwd=2, col=col["AxA"])
	points(0, 0, pch=16, col=col["add"], cex=3)
	lines(ellipse::ellipse(var(cbind(x, y))/nrow(dat), centre=c(mean(x), mean(y))), pch=1, lwd=3)
	legend("topleft", lty=c(0,0), pch=c(1, 16), col=c("black", col["add"]), legend=c("Data", "Additive"))
}




pdf("../results/Fig1BC.pdf", width=fig.width*2/3, height=fig.height, pointsize=fig.pointsize)
	layout(t(1:2))
	par(cex=1, mar=fig.mar) # layout changes the point size
	crossfig2(crosses.lines$Weight, main=weight.name)
	crossfig2(crosses.lines$Fitness, main=silique.name)
dev.off()


pdf("../results/FigS1.pdf", width=fig.width*2/3, height=fig.height, pointsize=fig.pointsize)
	layout(t(1:2))
	par(cex=1, mar=fig.mar) # layout changes the point size
	crossfig2(crosses.pops$Weight, main=paste0("Population: ", weight.name))
	crossfig2(crosses.pops$Fitness, main=paste0("Population: ", silique.name))
dev.off()

	dd.intrapop <- dd[dd$Mother_pop == dd$Father_pop,]
	dd.intrapop$cross <- with(dd.intrapop, paste0(Mother_pop, "-", sapply(strsplit(Mother_line, split="-"), function(x) x[2]), "x", sapply(strsplit(Father_line, split="-"), function(x) x[2])))

pdf("../results/FigS2a.pdf", width=fig.width, height=1.5*fig.height, pointsize=fig.pointsize)
	par(mar=fig.mar+c(2,0,0,0))
	bycross.Weight <- by(dd.intrapop$Weight, dd.intrapop$cross, FUN=c)
	bycross.pop <- sapply(strsplit(names(bycross.Weight), split="-"), FUN="[", 1)
	bycross.at  <- seq_along(bycross.pop) + cumsum(c(0,2*diff(as.numeric(factor(bycross.pop)))))
	boxplot(bycross.Weight, col=col.pops[bycross.pop], las=2, ylab=weight.name, at=bycross.at)
	legend("topright", lty=0, pch=22, pt.bg=col.pops, pt.cex=2, legend=names(col.pops), horiz=TRUE)
dev.off()

pdf("../results/FigS2b.pdf", width=fig.width, height=1.5*fig.height, pointsize=fig.pointsize)
	par(mar=fig.mar+c(2,0,0,0))
	bycross.Fitness <- by(dd.intrapop$Fitness, dd.intrapop$cross, FUN=c)
	bycross.pop <- sapply(strsplit(names(bycross.Fitness), split="-"), FUN="[", 1)
	bycross.at  <- seq_along(bycross.pop) + cumsum(c(0,2*diff(as.numeric(factor(bycross.pop)))))
	boxplot(bycross.Fitness, col=col.pops[bycross.pop], las=2, ylab=silique.name, at=bycross.at)
	#legend("topleft", lty=0, pch=15, col=col.pops, legend=names(col.pops), horiz=TRUE)
dev.off()

