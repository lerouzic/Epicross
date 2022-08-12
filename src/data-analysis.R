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
		P = sapply(lines, function(line) mean(dd[dd$Mother_line == line & dd$Father_line == line, trait])), 
		F1 = outer(lines, lines, function(line1, line2) mapply(line1, line2, FUN=function(l1, l2) mean(dd[dd$Mother_line == l1 & dd$Father_line == l2 & dd$Gen == "F1", trait]))),
		F2 = outer(lines, lines, function(line1, line2) mapply(line1, line2, FUN=function(l1, l2) mean(dd[dd$Mother_line == l1 & dd$Father_line == l2 & dd$Gen == "F2", trait])))
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


model.table <- do.call(rbind, list(
	A      = data.frame(A.pop=TRUE,  A.line=FALSE, D.pop=FALSE, D.line=FALSE, E.pop=FALSE, E.line=FALSE),
	lA     = data.frame(A.pop=FALSE, A.line=TRUE,  D.pop=FALSE, D.line=FALSE, E.pop=FALSE, E.line=FALSE),
	AlA    = data.frame(A.pop=TRUE,  A.line=TRUE,  D.pop=FALSE, D.line=FALSE, E.pop=FALSE, E.line=FALSE), 
	A.D    = data.frame(A.pop=TRUE,  A.line=FALSE, D.pop=TRUE , D.line=FALSE, E.pop=FALSE, E.line=FALSE),
	lA.D   = data.frame(A.pop=FALSE, A.line=TRUE,  D.pop=TRUE , D.line=FALSE, E.pop=FALSE, E.line=FALSE),
	lA.lD  = data.frame(A.pop=FALSE, A.line=TRUE,  D.pop=FALSE, D.line=TRUE,  E.pop=FALSE, E.line=FALSE),
	A.E    = data.frame(A.pop=TRUE,  A.line=FALSE, D.pop=FALSE, D.line=FALSE, E.pop=TRUE,  E.line=FALSE),
	lA.lE  = data.frame(A.pop=FALSE, A.line=TRUE,  D.pop=FALSE, D.line=FALSE, E.pop=FALSE, E.line=TRUE ),
	A.D.E  = data.frame(A.pop=TRUE,  A.line=FALSE, D.pop=TRUE , D.line=FALSE, E.pop=TRUE,  E.line=FALSE),
	lA.lD.lE=data.frame(A.pop=FALSE, A.line=TRUE,  D.pop=FALSE ,D.line=TRUE,  E.pop=FALSE, E.line=TRUE )
))

model.names <- setNames(nm=rownames(model.table))

#~ model.table <- model.table[1:2,]

#~ hglm.weight <- mclapply(model.names, function(mm) do.call(hglm.wrapper, c(as.list(model.table[mm,]), list(dd=dd, trait="Weight"))), mc.cores=mc.cores)
admb.weight <- mclapply(rownames(model.table), function(mm) do.call(admb.wrapper, c(as.list(model.table[mm,]), list(dd=dd, trait="Weight"))), mc.cores=mc.cores)

#~ hglm.fitness <- mclapply(model.names, function(mm) do.call(hglm.wrapper, c(as.list(model.table[mm,]), list(dd=dd, trait="Fitness"))), mc.cores=mc.cores)
admb.fitness <- mclapply(rownames(model.table), function(mm) do.call(admb.wrapper, c(as.list(model.table[mm,]), list(dd=dd, trait="Fitness"))), mc.cores=mc.cores)

#~ hglm.weight.summ <- do.call(rbind, lapply(hglm.weight, hglm.summ))
admb.weight.summ <- do.call(rbind, lapply(admb.weight, admb.summ))

#~ hglm.fitness.summ <- do.call(rbind, lapply(hglm.fitness, hglm.summ))
admb.fitness.summ <- do.call(rbind, lapply(admb.fitness, admb.summ))

bycross.weight <- t(apply(crosses.lines$Weight, 1, function(x) fit.2pop(x["P1"],x["P2"],x["F1"],x["F2"])))


library(viridis)

makeTransparent<-function(someColor, alpha=70)
{ # from https://stackoverflow.com/questions/8047668/transparent-equivalent-of-given-color
  newColor<-col2rgb(someColor)
  apply(newColor, 2, function(curcoldata){rgb(red=curcoldata[1], green=curcoldata[2],
    blue=curcoldata[3],alpha=alpha, maxColorValue=255)})
}

cols <- inferno(nrow(crosses.lines$Weight))
plot(NULL, xlim=c(0.25, 0.75), ylim=c(450, 1700), xlab="Weight", ylab="Fitness")
points(crosses.lines$Weight$P1[-24],  crosses.lines$Fitness$P1, col=cols, pch=4)
points(crosses.lines$Weight$P2[-24],  crosses.lines$Fitness$P2, col=cols, pch=4)
points(crosses.lines$Weight$F1[-24],  crosses.lines$Fitness$F1, col=cols, pch=1)
points(crosses.lines$Weight$F2[-24],  crosses.lines$Fitness$F2, col=cols, pch=16)

for (i in seq_along(cols)) {
	polygon(x=c(crosses.lines$Weight$P1[i], crosses.lines$Weight$F1[i], crosses.lines$Weight$P2[i], crosses.lines$Weight$F2[i], crosses.lines$Weight$P1[i]),
			y=c(crosses.lines$Fitness$P1[i], crosses.lines$Fitness$F1[i], crosses.lines$Fitness$P2[i], crosses.lines$Fitness$F2[i], crosses.lines$Fitness$P1[i]),
			border=cols[i], 
			col=makeTransparent(cols[i])
		)}



points((crosses.lines$Weight$P1 + crosses.lines$Weight$P2)/2, crosses.lines$Weight$F2, col=cols, pch=16)

abline(a=0, b=1, col="gray")


pdf("cross-fig.pdf", width=10, height=5)
	layout(t(1:2))
	plot(NULL, xlim=range(crosses.lines$Weight$sF1), ylim=range(crosses.lines$Weight$sF2),  xlab="relative F1", ylab="relative F2", main="weight", asp=1)
	points(crosses.lines$Weight$sF1, crosses.lines$Weight$sF2)
	abline(v=0, h=0, col="gray")
	abline(a=0, b=1, col="darkgreen", lwd=2)
	text(10,10,"No dominance", pos=3, srt=45, col="darkgreen")
	abline(a=0, b=0.5, col="darkviolet", lwd=2)
	text(10,5,"No epistasis", pos=3, srt=28, col="darkviolet")
	
	
	plot(NULL, xlim=range(crosses.lines$Fitness$sF1), ylim=range(crosses.lines$Fitness$sF2), xlab="relative F1", ylab="relative F2", main="fitness")
	points(crosses.lines$Fitness$sF1, crosses.lines$Fitness$sF2)
	abline(v=0, h=0, col="gray")
	abline(a=0, b=1, col="darkgreen", lwd=2)
	text(20,20,"No dominance", pos=3, srt=45, col="darkgreen")
	abline(a=0, b=0.5, col="darkviolet", lwd=2)
	text(20,10,"No epistasis", pos=3, srt=28, col="darkviolet")
dev.off()
