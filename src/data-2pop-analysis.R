source("../src/genet-model.R")

library(viridis)

extract.cross <- function(data, parent1, parent2) {
	list(
		P1 = data[(data$Mother_line == parent1 | data$Mother_pop == parent1) & (data$Father_line == parent1 | data$Father_pop == parent1),],
		P2 = data[(data$Mother_line == parent2 | data$Mother_pop == parent2) & (data$Father_line == parent2 | data$Father_pop == parent2),],
		F1 = data[(data$Mother_line == parent1 | data$Mother_pop == parent1 | data$Father_line == parent1 | data$Father_pop == parent1)
				  &(data$Mother_line == parent2 | data$Mother_pop == parent2 | data$Father_line == parent2 | data$Father_pop == parent2)
				  & data$Gen == "F1",],
		F2 = data[(data$Mother_line == parent1 | data$Mother_pop == parent1 | data$Father_line == parent1 | data$Father_pop == parent1)
				  &(data$Mother_line == parent2 | data$Mother_pop == parent2 | data$Father_line == parent2 | data$Father_pop == parent2)
				  & data$Gen == "F2",]
	)
}

model.selection.cross <- function(data, parent1, parent2, trait="Weight") {
	cc <- extract.cross(data, parent1, parent2)
	if (any(sapply(cc, nrow) == 0)) return(list(coef=NULL, AIC=NULL))
	
	model.a      <- fit.2pop.F2(cc$P1[,trait], cc$P2[,trait], cc$F1[,trait], cc$F2[,trait], what="a")
	model.a.d    <- fit.2pop.F2(cc$P1[,trait], cc$P2[,trait], cc$F1[,trait], cc$F2[,trait], what=c("a","d"))
	model.a.aa   <- fit.2pop.F2(cc$P1[,trait], cc$P2[,trait], cc$F1[,trait], cc$F2[,trait], what=c("a","aa"))
	model.a.d.aa <- fit.2pop.F2(cc$P1[,trait], cc$P2[,trait], cc$F1[,trait], cc$F2[,trait], what=c("a","d","aa"))
	
	list(
		coef = c(coef(model.a), coef(model.a.d), coef(model.a.aa), coef(model.a.d.aa)),
		AIC  = c(model.a=AIC(model.a), model.a.d=AIC(model.a.d), model.a.aa=AIC(model.a.aa), model.a.d.aa=AIC(model.a.d.aa))
	)
}

barplot.AkaikeW <- function(delta.AIC, what=c("Weight", "Fitness"), sort=FALSE) {
	layout(seq_along(what))
	par(mar=c(6, 4, 5, 1))
	for (ww in what) {
		AkW <- apply(delta.AIC[[ww]], 1, function(x) exp(-0.5*x)/sum(exp(-0.5*x)))
		if (sort)
			AkW <- AkW[,order(AkW["model.a",], decreasing=TRUE)]
		barplot(AkW, col=magma(4), las=2, ylab="Akaike Weight")
		title(ww, line=2.5)
		legend("top", horiz=TRUE, lwd=5, col=magma(4), legend=rownames(AkW), inset=c(0,-0.1,0), xpd=NA, bty="n", cex=0.8)
	}
}

dd <- read.table("../data/data_clean.txt", stringsAsFactors=FALSE)

#focus on pops
pops <- unique(c(dd$Mother_pop, dd$Father_pop))
cpops <- combn(pops, 2)
colnames(cpops) <- apply(cpops, 2, paste, collapse="x")
AICpop <- list(
	Weight  = t(mapply(cpops[1,], cpops[2,], FUN=function(p1, p2) { model.selection.cross(dd, p1, p2, trait="Weight")$AIC})), 
	Fitness = t(mapply(cpops[1,], cpops[2,], FUN=function(p1, p2) { model.selection.cross(dd, p1, p2, trait="Fitness")$AIC})))
delta.AICpop <- lapply(AICpop, function(x) t(apply(x, 1, function(xx) xx-min(xx))))

pdf("../results/F2AICpop.pdf", width=15, height=10)
	barplot.AkaikeW(delta.AICpop)
dev.off()

# focus on lines
lines <- unique(c(dd$Mother_line, dd$Father_line))
clines <- combn(lines, 2)
colnames(clines) <- apply(clines, 2, paste, collapse="x")
AICline <- list(
	Weight = do.call(rbind, mapply(clines[1,], clines[2,], FUN=function(p1, p2) { model.selection.cross(dd, p1, p2, trait="Weight")$AIC}, SIMPLIFY=FALSE)), 
	Fitness = do.call(rbind, mapply(clines[1,], clines[2,], FUN=function(p1, p2) { model.selection.cross(dd, p1, p2, trait="Fitness")$AIC}, SIMPLIFY=FALSE)))
delta.AICline <- lapply(AICline, function(x) t(apply(x, 1, function(xx) xx-min(xx))))

pdf("../results/F2AIClines.pdf", width=15, height=10)
	barplot.AkaikeW(delta.AICline)
dev.off()
