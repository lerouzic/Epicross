
# Sets the coefficient loadings (additive, dominance, 
# add x add, add x dom, and dom x dom)
# for parental, F1, and F2 populations. 

model1 <- rbind(
	P=c(a= 1, d=0,   aa=0,  ad=0,  dd=0),
	F1=c(a= 1, d=1,   aa=1,  ad=0,  dd=1),
	F2=c(a= 1, d=1/2, aa=1,  ad=0,  dd=1/4))
	
MODEL <- model1
