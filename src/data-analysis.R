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

draw_stars <- function(y1, y2, y3, p) {
	stars <- function(p) ifelse(p < 0.001, "***", ifelse(p < 0.01, "**", ifelse(p < 0.05, "*", "")))
	# Very specialized use; contrasts are expected to be in the right order
	ps <- stars(p)
	x1 <- c(-1, 2, 5, 8, 9, 4.5, 8, 8)
	x2 <- c(-1, 1, 4, 7, 8, 1.5, 1.5, 4.5)
	y  <- c(-1, y1, y1, y1, y1, y2, y3, y2)
	for (i in seq_along(p)) {
		if (i > 1 && ps[i] != "") {
			arrows(x0=x1[i]-0.1, x1=x2[i]+0.1, y0=y[i], angle=90, length=0.005, code="3")
			text(0.5*x1[i]+0.5*x2[i], y[i], ps[i], adj=c(0.5,0.2))
		}
	}
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
	crossfig2(crosses.pops$Weight, main=paste0("Pop: ", weight.name))
	crossfig2(crosses.pops$Fitness, main=paste0("Pop: ", silique.name))
dev.off()

	dd.intrapop <- dd[dd$Mother_pop == dd$Father_pop,]
	dd.intrapop$cross <- with(dd.intrapop, paste0(Mother_pop, "-", sapply(strsplit(Mother_line, split="-"), function(x) x[2]), "x", sapply(strsplit(Father_line, split="-"), function(x) x[2])))

linesort <- function(linenames) {
	ss <- strsplit(linenames, split="[-x]")
	sss <- split(ss, sapply(ss, "[", 1))
	o <- unlist(lapply(sss, function(x) x[order(sapply(x, "[", 2) == sapply(x, "[", 3))]), recursive=FALSE)
	paste0(sapply(o, "[", 1), "-", sapply(o, "[", 2), "x", sapply(o, "[", 3))
}

pdf("../results/FigS2a.pdf", width=fig.width, height=1.5*fig.height, pointsize=fig.pointsize)
	par(mar=fig.mar+c(2,0,0,0))
	bycross.Weight <- by(dd.intrapop$Weight, dd.intrapop$cross, FUN=c)
	bycross.Weight <-  bycross.Weight[linesort(names(bycross.Weight))]
	maternal.pop <- sapply(strsplit(names(bycross.Weight), split="[-x]"), function(x) length(unique(x))) == 2
	bycross.pop <- sapply(strsplit(names(bycross.Weight), split="-"), FUN="[", 1)
	bycross.at  <- seq_along(bycross.pop) + cumsum(c(0,2*diff(as.numeric(factor(bycross.pop)))))
	boxplot(bycross.Weight, col=col.pops[bycross.pop], las=2, ylab=weight.name, at=bycross.at, border=ifelse(maternal.pop, "black", "grey40"))
	legend("topright", lty=0, pch=22, pt.bg=col.pops, pt.cex=2, legend=names(col.pops), horiz=TRUE)
dev.off()

pdf("../results/FigS2b.pdf", width=fig.width, height=1.5*fig.height, pointsize=fig.pointsize)
	par(mar=fig.mar+c(2,0,0,0))
	bycross.Fitness <- by(dd.intrapop$Fitness, dd.intrapop$cross, FUN=c)
	bycross.Fitness <-  bycross.Fitness[linesort(names(bycross.Fitness))]
	maternal.pop <- sapply(strsplit(names(bycross.Fitness), split="[-x]"), function(x) length(unique(x))) == 2
	bycross.pop <- sapply(strsplit(names(bycross.Fitness), split="-"), FUN="[", 1)
	bycross.at  <- seq_along(bycross.pop) + cumsum(c(0,2*diff(as.numeric(factor(bycross.pop)))))
	boxplot(bycross.Fitness, col=col.pops[bycross.pop], las=2, ylab=silique.name, at=bycross.at, border=ifelse(maternal.pop, "black", "grey40"))
	#legend("topleft", lty=0, pch=15, col=col.pops, legend=names(col.pops), horiz=TRUE)
dev.off()


require(multcomp)

pdf("../results/FigS4.pdf", width=fig.width, height=fig.height, pointsize=fig.pointsize)
	layout(t(1:2))
	par(cex=1, mar=fig.mar+c(2,0,0,0)) 
	
	mother_eco <- ifelse(dd$Mother_pop %in% names(col.pops)[1:3], "High", "Low")
	father_eco <- ifelse(dd$Father_pop %in% names(col.pops)[1:3], "High", "Low")
	dd$Type <- factor(
		ifelse(dd$Mother_pop == dd$Father_pop, 
		       ifelse(dd$Mother_line == dd$Father_line, paste0("Self_", mother_eco), paste0("Within_", mother_eco)),
		       ifelse(mother_eco == father_eco, paste0(mother_eco, " x ", father_eco), "High x Low")),
		levels = c("Self_High", "Self_Low", "Within_High", "Within_Low", "High x High", "High x Low", "Low x Low"))
	
	# Contrasts
	Type.ct <- rbind(
		`Intercept`     = c(1, 0, 0, 0, 0, 0, 0),
		`Self_Low-High` = c(-1,1, 0, 0, 0, 0, 0),
		`Within_Low-High`=c(0, 0,-1, 1, 0, 0, 0),
		`HL-HH`          =c(0, 0, 0, 0,-1, 1, 0),
		`LL-HL`          =c(0, 0, 0, 0, 0,-1, 1),
		`Within-Self`    =c(-0.5,-0.5,0.5,0.5,0,0,0),
		`Between-Self`   =c(-0.5,-0.5,0,0,1/3,1/3,1/3),
		`Between-Within` =c(0, 0,-0.5,-0.5,1/3,1/3,1/3))
	
	llw <- lm(Weight ~ 0 + Type, data=dd)
	llf <- lm(Fitness ~ 0 + Type, data=dd)
	
	pvalw <- summary(glht(llw, linfct=mcp(Type=Type.ct)), test=adjusted("none"))$test$pvalues
	pvalf <- summary(glht(llf, linfct=mcp(Type=Type.ct)), test=adjusted("none"))$test$pvalues

	boxplot(dd$Weight ~ dd$Type, at=c(1:2, 4:5, 7:9), ylim=c(0, 1.5), col=col.pops[c(2,5,2,5,2,NA,5)], xlab="", ylab=weight.name, las=2)
	draw_stars(y1=1.3, y2=1.45, y3=1.5, p=pvalw)
	legend("topleft", fill=col.pops[c(2,5,NA)], legend=c("High", "Low", "Hybrid"), bty="n", horiz=TRUE)
	
	boxplot(dd$Fitness ~ dd$Type, at=c(1:2, 4:5, 7:9), ylim=c(0, 3100), col=col.pops[c(2,5,2,5,2,NA,5)], xlab="", ylab=silique.name, las=2)
	draw_stars(y1=2600, y2=2800, y3=3000, p=pvalf)

dev.off()

