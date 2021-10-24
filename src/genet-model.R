
# Sets the coefficient loadings (additive, dominance, 
# add x add, add x dom, and dom x dom)
# for parental, F1, and F2 populations. 

model1 <- rbind(
	P1=c(a1=1,   a2=0,   d=0,   aa=0,    ad=0, dd=0),
	P2=c(a1=0,   a2=1,   d=0,   aa=0,    ad=0, dd=0), 
	F1=c(a1=1/2, a2=1/2, d=1,   aa=-1/2, ad=0, dd=1/2),
	F2=c(a1=1/2, a2=1/2, d=1/2, aa=-1/2, ad=0, dd=1/8))
	
MODEL <- model1
