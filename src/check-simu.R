source("../src/data-simulator.R")
source("../src/admb-helper.R")

library(parallel)
mc.cores <- if (Sys.info()["sysname"] == "Windows") 1 else min(8, detectCores()-1)

.getparfit <- function(fit, dim1, dim2, param) {
	outer(setNames(nm=dim1), setNames(nm=dim2), function(d1, d2) mapply(d1, d2, FUN=function(dd1, dd2) fit[[dd1]][[dd2]][[param]]))
}

# Population effects
ss.pop.param <- list(
	pA   = list(varA.pop=1),
	pAD  = list(varA.pop=1, varD.pop=1),
	pADd = list(varA.pop=1, varD.pop=1, meanD.pop=-1),
	pAE  = list(varA.pop=1, varE.pop=1),
	pAEe = list(varA.pop=1, varE.pop=1, meanE.pop=-1),
	pADE = list(varA.pop=1, varD.pop=1, varE.pop=1)
)

mm.pop.param <- list(
	pA     = list(A.pop=TRUE,  A.line=FALSE, D.pop=FALSE, D.line=FALSE, E.pop=FALSE, E.line=FALSE),
	pAD    = list(A.pop=TRUE,  A.line=FALSE, D.pop=TRUE,  D.line=FALSE, E.pop=FALSE, E.line=FALSE),
	pAE    = list(A.pop=TRUE,  A.line=FALSE, D.pop=FALSE, D.line=FALSE, E.pop=TRUE,  E.line=FALSE),
	pADE   = list(A.pop=TRUE,  A.line=FALSE, D.pop=TRUE,  D.line=FALSE, E.pop=TRUE,  E.line=FALSE))


fit.pop <- lapply(ss.pop.param, function(ss) { # do not parallelize this one (running the same admb model in parallel is a bad idea)
	mclapply(mm.pop.param, function(mm) {
		admb.summ(do.call(admb.wrapper, c(list(dd=do.call(simulate.data, ss), trait="Phen"), mm)))
	}, mc.cores=mc.cores)
})

par.pop <- list(
	Vpop.a = .getparfit(fit.pop, names(mm.pop.param), names(ss.pop.param), "Vpop.a"), 
	Vpop.d = .getparfit(fit.pop, names(mm.pop.param), names(ss.pop.param), "Vpop.d"), 
	Vpop.e = .getparfit(fit.pop, names(mm.pop.param), names(ss.pop.param), "Vpop.e")
) 
	
