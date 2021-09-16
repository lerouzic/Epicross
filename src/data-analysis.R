library(parallel)

mc.cores <- min(8, detectCores() - 1)
#~ mc.cores <- 1

source("../src/admb-helper.R")
source("../src/hglm-helper.R")

dd <- read.table("../data/data_clean.txt", stringsAsFactors=FALSE)

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

hglm.weight <- mclapply(model.names, function(mm) do.call(hglm.wrapper, c(as.list(model.table[mm,]), list(dd=dd, trait="Weight"))), mc.cores=mc.cores)
admb.weight <- mclapply(model.names, function(mm) do.call(admb.wrapper, c(as.list(model.table[mm,]), list(dd=dd, trait="Weight"))), mc.cores=mc.cores)

hglm.fitness <- mclapply(model.names, function(mm) do.call(hglm.wrapper, c(as.list(model.table[mm,]), list(dd=dd, trait="Fitness"))), mc.cores=mc.cores)
admb.fitness <- mclapply(model.names, function(mm) do.call(admb.wrapper, c(as.list(model.table[mm,]), list(dd=dd, trait="Fitness"))), mc.cores=mc.cores)

hglm.weight.summ <- do.call(rbind, lapply(hglm.weight, hglm.summ))
admb.weight.summ <- do.call(rbind, lapply(admb.weight, admb.summ))

hglm.fitness.summ <- do.call(rbind, lapply(hglm.fitness, hglm.summ))
admb.fitness.summ <- do.call(rbind, lapply(admb.fitness, admb.summ))

