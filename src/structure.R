
model1 <- rbind(
	P=c(a= 1, d=0,   aa=0,  ad=0,  dd=0),
	F1=c(a= 1, d=1,   aa=1,  ad=0,  dd=1),
	F2=c(a= 1, d=1/2, aa=1,  ad=0,  dd=1/4))
	
MODEL <- model1

.check.par <- function(pop.par1, pop.par2, gen, effect) {
	stopifnot(
		gen %in% c("F1","F2"),
		is.vector(pop.par1), is.vector(pop.par2),
		length(pop.par1) == length(pop.par2),
		effect %in% colnames(MODEL))
}

structure.matrix.1D <- function(pop.par1, pop.par2, gen, effect="a") {
	.check.par(pop.par1, pop.par2, gen, effect)
	
	lv <- gtools::mixedsort(levels(factor(c(as.character(pop.par1), as.character(pop.par2)))))
	
	X <- matrix(0, nrow=length(pop.par1), ncol=length(lv))
	colnames(X) <- lv
	gen <- ifelse(pop.par1 == pop.par2, "P", gen)
	
	for (i in seq_along(pop.par1))
		X[i, as.character(pop.par1[i])] <- X[i, as.character(pop.par1[i])] + MODEL[gen[i], effect]
	for (i in seq_along(pop.par2))
		X[i, as.character(pop.par2[i])] <- X[i, as.character(pop.par2[i])] + MODEL[gen[i], effect]
	# Removing empty columns
	X[,colSums(abs(X))>0] 
}

structure.matrix.2D <- function(pop.par1, pop.par2, gen, effect="d") {
	.check.par(pop.par1, pop.par2, gen, effect)
	
	# Not very efficient algorithm. 
	lv <- gtools::mixedsort(levels(factor(c(as.character(pop.par1), as.character(pop.par2)))))
	lv.inter <- outer(lv, lv, paste, sep="-")
	lv.inter <- c(lv.inter[upper.tri(lv.inter, diag=TRUE)])
	
	X <- matrix(0, nrow=length(pop.par1), ncol=length(lv.inter))
	colnames(X) <- lv.inter
	gen <- ifelse(pop.par1 == pop.par2, "P", gen)
	
	for (i in seq_along(pop.par1)) {
		nm <- paste(as.character(pop.par1[i]), as.character(pop.par2[i]), sep="-")
		if (!nm %in% colnames(X))
			nm <- paste(as.character(pop.par2[i]), as.character(pop.par1[i]), sep="-")
		X[i, nm] <- X[i, nm] + MODEL[gen[i], effect]
	}
	# Removing empty columns
	X[,colSums(abs(X))>0]
}

additive.effects <- function(pop.par1, pop.par2, gen) {
	list(
		X=NULL,
		Z=structure.matrix.1D(pop.par1, pop.par2, gen, "a")
	)
}

dominance.effects <- function(pop.par1, pop.par2, gen) {
	Zdom <- structure.matrix.2D(pop.par1, pop.par2, gen, "d")
	list(
		X=as.matrix(rowSums(Zdom)),
		Z=Zdom
	)
}

AxA.effects <- function(pop.par1, pop.par2, gen) {
	Zaa <- structure.matrix.2D(pop.par1, pop.par2, gen, "aa")
	list(
		X=as.matrix(rowSums(Zaa)),
		Z=Zaa
	)
}

DxD.effects <- function(pop.par1, pop.par2, gen) {
	Zdd <- structure.matrix.2D(pop.par1, pop.par2, gen, "dd")
	list(
		X=as.matrix(rowSums(Zdd)),
		Z=Zdd
	)
}
