source("../src/data-simulator.R")
source("../src/admb-helper.R")

library(parallel)
library(R.utils)

num.tests <- 21

mc.cores <- if (Sys.info()["sysname"] == "Windows") 1 else min(num.tests, detectCores()-1)

cache.qq.pop <- "../cache/fit-qq-pop.rds"
cache.qq.lin <- "../cache/fit-qq-lin.rds"


par.translate <- c(
	varA.pop = "Vpop.a",
	varD.pop = "Vpop.d",
	varE.pop = "Vpop.e",
	meanD.pop= "pop.d",
	meanE.pop= "pop.e",
	varA.line= "Vlin.a",
	varD.line= "Vlin.d",
	varE.line= "Vlin.e",
	meanD.line="lin.d",
	meanE.line="lin.e")

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

mm.lin.param <- list(
	pADE   = list(A.pop=FALSE, A.line=TRUE,  D.pop=FALSE, D.line=TRUE,  E.pop=FALSE, E.line=TRUE))

#~ fit.pop <- lapply(ss.pop.param, function(ss) { # do not parallelize this one (running the same admb model in parallel is a bad idea)
#~ 	mclapply(mm.pop.param, function(mm) {
#~ 		admb.summ(do.call(admb.wrapper, c(list(dd=do.call(simulate.data, ss), trait="Phen"), mm)))
#~ 	}, mc.cores=mc.cores)
#~ })

#~ par.pop <- list(
#~ 	Vpop.a = .getparfit(fit.pop, names(mm.pop.param), names(ss.pop.param), "Vpop.a"), 
#~ 	Vpop.d = .getparfit(fit.pop, names(mm.pop.param), names(ss.pop.param), "Vpop.d"), 
#~ 	Vpop.e = .getparfit(fit.pop, names(mm.pop.param), names(ss.pop.param), "Vpop.e")
#~ ) 
	

qq.pop.param <- list(varA.pop=1, varD.pop=1, meanD.pop=0, varE.pop=1, meanE.pop=0)

if (!file.exists(cache.qq.pop)){
	fit.qq.pop <- lapply(setNames(nm=names(qq.pop.param)), function(pp) {
		xx <- seq(-1, 1, length.out=num.tests)
		if (grepl("var", pp)) xx <- xx+1 # No negative numbers for variances 
		mclapply(setNames(xx, nm=as.character(round(xx, digits=3))), function(ppp) {
			my.qq <- qq.pop.param
			my.qq[[pp]] <- ppp
			withTimeout(
				admb.summ(do.call(admb.wrapper, c(list(dd=do.call(simulate.data, my.qq), trait="Phen"), mm.pop.param[["pADE"]]))),
				timeout=100, onTimeout="warning")
		}, mc.cores=mc.cores)
	})
	saveRDS(fit.qq.pop, file=cache.qq.pop, version=2)
}
fit.qq.pop <- readRDS(cache.qq.pop)


qq.lin.param <- list(varA.line=1, varD.line=1, meanD.line=0, varE.line=1, meanE.line=0)

if (!file.exists(cache.qq.lin)){
	fit.qq.lin <- lapply(setNames(nm=names(qq.lin.param)), function(pp) {
		xx <- seq(-1, 1, length.out=num.tests)
		if (grepl("var", pp)) xx <- xx+1 # No negative numbers for variances 
		mclapply(setNames(xx, nm=as.character(round(xx, digits=3))), function(ppp) {
			my.qq <- qq.lin.param
			my.qq[[pp]] <- ppp
			withTimeout(
				admb.summ(do.call(admb.wrapper, c(list(dd=do.call(simulate.data, my.qq), trait="Phen"), mm.lin.param[["pADE"]]))),
				timeout=1000, onTimeout="warning")
		}, mc.cores=mc.cores)
	})
	saveRDS(fit.qq.lin, file=cache.qq.lin, version=2)
}
fit.qq.lin <- readRDS(cache.qq.lin)





layout(1:length(fit.qq.pop))

for (nn in names(fit.qq.pop)) {
	xx <- as.numeric(names(fit.qq.pop[[nn]]))
	yy <- sapply(fit.qq.pop[[nn]], function(ll) if(is.list(ll)) return(ll[[par.translate[nn]]]) else NA)
	plot(xx, yy, xlab=paste0("Simulated ", nn), ylab=paste0("Estimated ", nn))
	abline(a=0, b=1, col="gray")
}

layout(1:length(fit.qq.lin))

for (nn in names(fit.qq.lin)) {
	xx <- as.numeric(names(fit.qq.lin[[nn]]))
	yy <- sapply(fit.qq.lin[[nn]], function(ll) if(is.list(ll)) return(ll[[par.translate[nn]]]) else NA)
	plot(xx, yy, ylim=range(yy, na.rm=TRUE), xlab=paste0("Simulated ", nn), ylab=paste0("Estimated ", nn))
	abline(a=0, b=1, col="gray")
}
