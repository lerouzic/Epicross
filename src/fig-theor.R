source("../src/genet-model.R")

get.pop.raw <- function(pop, a1, a2, d, aa, ad, dd) {
	a1*MODEL[pop, "a1"] +a2*MODEL[pop, "a2"] + d*MODEL[pop, "d"] + aa*MODEL[pop, "aa"] + ad*MODEL[pop, "ad"] + dd*MODEL[pop, "dd"]
}

get.pop <- function(pop, a1, a2, d, e) {
	get.pop.raw(pop, a1, a2, d, a1*a2*e, (a1+a2)*d*e, d*d*e)
}

get.d.e <- function(P1=1, P2=-1, F1, F2) {
	F1 <- (F1-(0.5*P1+0.5*P2))/(P1-P2)*2
	F2 <- (F2-(0.5*P1+0.5*P2))/(P1-P2)*2
	optim(c(d=0, e=0), fn=function(pp) (pp["d"]-(1/2)*pp["e"] + (1/2)*pp["e"]*pp["d"]^2-F1)^2 + (pp["d"]/2-pp["e"]/2+(1/8)*pp["e"]*pp["d"]^2-F2)^2)$par
}

layout(1)

plot(NULL, xlim=c(-1, 1), ylim=c(-1, 1), ylab="Phenotype", xlab="Dominance")

abline(h=c(-1,1), col="black")

curve(get.pop("F1", a1=1, a2=-1, d=x, e=0), col="green", add=TRUE)
curve(get.pop("F2", a1=1, a2=-1, d=x, e=0), col="green", lty=2, add=TRUE)

curve(get.pop("F1", a1=1, a2=-1, d=x, e=1), col="blue", add=TRUE)
curve(get.pop("F2", a1=1, a2=-1, d=x, e=1), col="blue", lty=2, add=TRUE)

curve(get.pop("F1", a1=1, a2=-1, d=x, e=-1), col="red", add=TRUE)
curve(get.pop("F2", a1=1, a2=-1, d=x, e=-1), col="red", lty=2, add=TRUE)

legend("topleft", lty=1, col=c("blue","green","red"), legend=c("e=1", "e=0", "e=-1"), bg="white")
legend("bottomright", lty=c(1,2), legend=c("F1","F2"), bg="white")



ff1 <- seq(-0.8, 0.8, length.out=41)
ff2 <- seq(-0.8, 0.8, length.out=41)

dd <- outer(ff1, ff2, function(ff1, ff2) mapply(ff1, ff2, FUN=function(fff1, fff2) get.d.e(F1=fff1, F2=fff2)["d"]))
ee <- outer(ff1, ff2, function(ff1, ff2) mapply(ff1, ff2, FUN=function(fff1, fff2) get.d.e(F1=fff1, F2=fff2)["e"]))

layout(t(1:2))
contour(ff1, ff2, dd, xlab="F1", ylab="F2", main="Dominance")
contour(ff1, ff2, ee, xlab="F1", ylab="F2", main="Epistasis")


# Multilocus approach

.which.a  <- function(l) which(sapply(strsplit(noia:::effectsNamesGeneral(l), split=""), function(x) sum(x=="a") == 1 & sum(x=="d") == 0))
.which.d  <- function(l) which(sapply(strsplit(noia:::effectsNamesGeneral(l), split=""), function(x) sum(x=="a") == 0 & sum(x=="d") == 1))
.which.aa <- function(l) which(sapply(strsplit(noia:::effectsNamesGeneral(l), split=""), function(x) sum(x=="a") == 2 & sum(x=="d") == 0))
.which.ad <- function(l) which(sapply(strsplit(noia:::effectsNamesGeneral(l), split=""), function(x) sum(x=="a") == 1 & sum(x=="d") == 1))
.which.dd <- function(l) which(sapply(strsplit(noia:::effectsNamesGeneral(l), split=""), function(x) sum(x=="a") == 0 & sum(x=="d") == 2))

.F2frq <- function(l) { ans <- 1; for (i in 1:l) ans <- ans %x% c(1/4, 1/2, 1/4); setNames(c(ans), nm=noia:::genNames(l)) }

.makeS <- function(ref, l) { ans <- 1; for (i in 1:l) ans <- ans %x% noia:::Sloc(ref); rownames(ans) <- c(noia:::genNames(l)); colnames(ans) <- noia:::effectsNamesGeneral(l); ans }

.genetmod.pop <- function(ref, l, pop) { 
	S <- .makeS(ref,l)
	vv <- setNames(rep(0, 3^l), nm=noia:::genNames(l))
	if (pop == "P1")
			vv[paste(rep('1',l),collapse="")] <- 1
		else if (pop == "P2")
			vv[paste(rep('3',l),collapse="")] <- 1
		else if (pop == "F1")
			vv[paste(rep('2',l),collapse="")] <- 1
		else if (pop == "F2")
			vv <- .F2frq(l)
	c(a=sum(S[,.which.a(l)] * vv), d=sum(S[,.which.d(l)] * vv), aa=sum(S[,.which.aa(l)] * vv), ad=sum(S[,.which.ad(l)] * vv), dd=sum(S[,.which.dd(l)] * vv))
}

.genetmod <- function(ref, l) { rbind(P1=.genetmod.pop(ref, l, "P1"), P2=.genetmod.pop(ref, l, "P2"), F1=.genetmod.pop(ref, l, "F1"), F2=.genetmod.pop(ref, l ,"F2")) }
