
source("../src/genet-model.R")

simulate.data <- function(N.pop=6, N.lines.per.pop=4, N.cross.within=1, N.cross.out=1, N.F1=10, N.F2=20, 
                          mu=0, varA.pop=0, varA.line=0, 
                          meanD.pop=0, meanD.line=0, varD.pop=0, varD.line=0, 
                          meanE.pop=0, meanE.line=0, varE.pop=0, varE.line=0, 
                          varR=1) {
	.wl <- function(cc) which(letters == cc)
	.wL <- function(cc) which(LETTERS == cc)
							  
	stopifnot(	N.pop > 1, N.lines.per.pop > 2, 
				N.cross.within > 0, N.cross.within < N.lines.per.pop, 
				N.cross.out > 0, N.cross.out < min(N.pop, N.lines.per.pop))
	
	# Names of random effect vectors
	pop.names        <- LETTERS[1:N.pop]
	line.names       <- c(t(outer(pop.names, letters[1:N.lines.per.pop], paste, sep="")))
	inter.pop.names  <- outer(pop.names, pop.names, paste, sep="")
	inter.pop.names  <- c(inter.pop.names[upper.tri(inter.pop.names, diag=TRUE)])
	inter.line.names <- outer(line.names, line.names, paste, sep="")
	inter.line.names <- c(inter.line.names[upper.tri(inter.line.names, diag=TRUE)])

	# Values of random effect vectors
	pops    <- setNames(rnorm(N.pop,                 0,             sqrt(varA.pop)),   nm=pop.names)
	lines   <- setNames(rnorm(N.pop*N.lines.per.pop, 0 ,            sqrt(varA.line)),  nm=line.names)
	dom.pops<- setNames(rnorm(length(inter.pop.names), meanD.pop,   sqrt(varD.pop)),   nm=inter.pop.names)
	dom     <- setNames(rnorm(length(inter.line.names), meanD.line, sqrt(varD.line)),  nm=inter.line.names)
	eps.pops<- setNames(rnorm(length(inter.pop.names), meanE.pop,   sqrt(varE.pop)),   nm=inter.pop.names)
	eps     <- setNames(rnorm(length(inter.line.names), meanE.line, sqrt(varE.line)),  nm=inter.line.names)
        # No interaction with itself
	dom.pops[paste(pop.names, pop.names, sep="")] <- 0
	dom     [paste(line.names,line.names,sep="")] <- 0
	eps.pops[paste(pop.names, pop.names, sep="")] <- 0
	eps     [paste(line.names,line.names,sep="")] <- 0
	
	# Self crosses
	moth.self <- fath.self <- line.names
	
	# Within-pop crosses
	moth.within <- rep(line.names, each=N.cross.within)
	fath.within <- unlist(sapply(strsplit(line.names, split=""), function(ll) {
						let <- letters[(.wl(ll[2]) + 1:(N.cross.within)) %% N.lines.per.pop]
						if (length(let) == 0) let <- letters[N.lines.per.pop]
						paste0(ll[1], let) }))
	
	# Among-pop crosses
	moth.among <- rep(line.names, each=N.cross.out)
	fath.among <- unlist(sapply(strsplit(line.names, split=""), function(ll) {
						letP <- LETTERS[(.wL(ll[1]) + 1:(N.cross.out)) %% N.pop]
						if (length(letP) == 0) letP <- LETTERS[N.pop]
						let <- letters[(.wl(ll[2]) + 0:(N.cross.out-1)) %% N.lines.per.pop]
						if (length(let) == 0) let <- letters[N.lines.per.pop]
						paste0(letP, let)}))
	
	moth <- c(rep(moth.self, each=N.F1), rep(moth.within, each=N.F1+N.F2), rep(moth.among, each=N.F1+N.F2))
	fath <- c(rep(fath.self, each=N.F1), rep(fath.within, each=N.F1+N.F2), rep(fath.among, each=N.F1+N.F2))
	inter<- paste(ifelse(moth < fath, moth, fath), ifelse(moth < fath, fath, moth), sep="")
	
	
	moth.pop <- sapply(strsplit(moth, split=""), "[", 1)
	fath.pop <- sapply(strsplit(fath, split=""), "[", 1)
	inter.pop<- paste(ifelse(moth.pop < fath.pop, moth.pop, fath.pop), ifelse(moth.pop < fath.pop, fath.pop, moth.pop), sep="")

	gen      <- c(	rep("F1", N.pop*N.lines.per.pop*N.F1),
					rep("F1", N.pop*N.lines.per.pop*N.cross.within*N.F1), rep("F2", N.pop*N.lines.per.pop*N.cross.within*N.F2),
					rep("F1", N.pop*N.lines.per.pop*N.cross.out*N.F1), rep("F2", N.pop*N.lines.per.pop*N.cross.out*N.F2))
					
	# Model coefficients from genet-model.R
	coef.a      <- ifelse(moth == fath,         MODEL["P","a"], ifelse(gen == "F1", MODEL["F1","a"], MODEL["F2","a"]))
	coef.a.pop  <- ifelse(moth.pop == fath.pop, MODEL["P","a"], ifelse(gen == "F1", MODEL["F1","a"], MODEL["F2","a"]))
	coef.d      <- ifelse(moth == fath,         MODEL["P","d"], ifelse(gen == "F1", MODEL["F1","d"], MODEL["F2","d"]))
	coef.d.pop  <- ifelse(moth.pop == fath.pop, MODEL["P","d"], ifelse(gen == "F1", MODEL["F1","d"], MODEL["F2","d"]))
	coef.aa     <- ifelse(moth == fath,         MODEL["P","aa"],ifelse(gen == "F1", MODEL["F1","aa"],MODEL["F2","aa"]))
	coef.aa.pop <- ifelse(moth.pop == fath.pop, MODEL["P","aa"],ifelse(gen == "F1", MODEL["F1","aa"],MODEL["F2","aa"]))
	coef.ad     <- ifelse(moth == fath,         MODEL["P","ad"],ifelse(gen == "F1", MODEL["F1","ad"],MODEL["F2","ad"]))
	coef.ad.pop <- ifelse(moth.pop == fath.pop, MODEL["P","ad"],ifelse(gen == "F1", MODEL["F1","ad"],MODEL["F2","ad"]))
	coef.dd     <- ifelse(moth == fath,         MODEL["P","dd"],ifelse(gen == "F1", MODEL["F1","dd"],MODEL["F2","dd"]))
	coef.dd.pop <- ifelse(moth.pop == fath.pop, MODEL["P","dd"],ifelse(gen == "F1", MODEL["F1","dd"],MODEL["F2","dd"]))
					
	# Predicted phenotype
	pred  <- mu +
	         coef.a.pop*(pops[moth.pop] + pops[fath.pop]) +
	         coef.a    *(lines[moth]    + lines[fath])    +
	         coef.d.pop*dom.pops[inter.pop] + 
	         coef.d    *dom[inter]          +
	         coef.aa.pop*pops[moth.pop]*pops[fath.pop]*eps.pops[inter.pop] +
	         coef.aa    *lines[moth]   *lines[fath]   *eps[inter]     +
	         coef.ad.pop*(pops[moth.pop] + pops[fath.pop])*dom.pops[inter.pop]*eps.pops[inter.pop] +
	         coef.ad    *(lines[moth]    + lines[fath]   )*dom[inter]         *eps[inter]     +
	         coef.dd.pop*(dom.pops[inter.pop]^2)*eps.pops[inter.pop] +
	         coef.dd    *(dom[inter]^2)         *eps[inter]

	# Simulated dataset
	ans <- data.frame(
		Mother_line = moth, 
		Mother_pop  = moth.pop,
		Father_line = fath,
		Father_pop  = fath.pop, 
		Gen         = gen, 
		Phen        = rnorm(length(pred), pred, sqrt(varR))
	)
}
