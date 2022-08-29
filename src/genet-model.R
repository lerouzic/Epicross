
# Sets the coefficient loadings (additive, dominance, 
# add x add, add x dom, and dom x dom)
# for parental, F1, and F2 populations. 

model1 <- rbind(
	P1=c(a1=1,   a2=0,   d=0,   aa=0,    ad=0, dd=0),
	P2=c(a1=0,   a2=1,   d=0,   aa=0,    ad=0, dd=0), 
	F1=c(a1=1/2, a2=1/2, d=1,   aa=-1/2, ad=0, dd=1/2),
	F2=c(a1=1/2, a2=1/2, d=1/2, aa=-1/2, ad=0, dd=1/8))

model2 <- rbind(
	P=c(a=1,    d=0,   aa=0,    ad=0, dd=0),
	F1=c(a=1/2, d=1,   aa=-1/2, ad=0, dd=1/2),
	F2=c(a=1/2, d=1/2, aa=-1/2, ad=0, dd=1/8))
	
MODEL <- model2



fit.2pop <- function(P1, P2, F1, F2) {
	ans <- try(optim(par=c(zR=unname((P1+P2))/2, a=unname((P2-P1))/2, d=unname(F1-(P1+P2)/2), e=0), 
	      fn=function(pp) 
	      	(P1 - (pp["zR"] - MODEL["P","a"] * pp["a"])) ^ 2 +
	      	(P2 - (pp["zR"] + MODEL["P","a"] * pp["a"])) ^ 2 +
	      	(F1 - (pp["zR"] + MODEL["F1","d"] * pp["d"] + MODEL["F1","dd"] * pp["e"] * pp["d"]^2  + MODEL["F1","aa"] * pp["e"] * pp["a"]^2 )) ^ 2 +
	      	(F2 - (pp["zR"] + MODEL["F2","d"] * pp["d"] + MODEL["F2","dd"] * pp["e"] * pp["d"]^2  + MODEL["F2","aa"] * pp["e"] * pp["a"]^2 )) ^ 2
	      ))
	if (class(ans) == "try-error") 
		return(rep(NA, 4))
	else
		ans$par
}

fit.2pop.noDD <- function(P1, P2, F1, F2) {
	ans <- c(
		zR = (P2+P1)/2,
		a  = (P2-P1)/2)
	F1s <- (F1-ans["zR"])/ans["a"]
	F2s <- (F2-ans["zR"])/ans["a"]
	ds <- 2*(F1s-F2s)
	es <- 2*(F1s-2*F2s)
	ans["d"] <- ds*ans["a"]
	ans["e"] <- es/ans["a"]
	ans
}
