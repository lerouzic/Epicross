library(hglm)
source("../src/hglm-copy.R")

library(R2admb)

dd <- read.table("../data/data_clean.txt", stringsAsFactors=FALSE)

source("../src/structure.R")

hglm.wrapper <- function(dd, trait="Weight", A.pop=TRUE, A.line=FALSE, D.pop=TRUE, D.line=FALSE, AA.pop=TRUE, AA.line=FALSE, DD.pop=TRUE, DD.line=FALSE) {
	
	X <- matrix(rep(1, nrow(dd)))
	
	listX <- list(mu=X)
	listZ <- list()
		
	if (A.pop) {
		ssA <- additive.effects(dd$Mother_pop, dd$Father_pop, dd$Gen)
		listZ <- c(listZ, list(pop.a=ssA$Z))
	}
	
	if (A.line) {
		ssA <- additive.effects(dd$Mother_line, dd$Father_line, dd$Gen)
		listZ <- c(listZ, list(lin.a=ssA$Z))
	}
	
	if (D.pop) {
		ssD <- dominance.effects(dd$Mother_pop, dd$Father_pop, dd$Gen)
		listX <- c(listX, list(pop.d=ssD$X))
		listZ <- c(listZ, list(pop.d=ssD$Z))
	}
	
	if (D.line) {
		ssD <- dominance.effects(dd$Mother_line, dd$Father_line, dd$Gen)
		listX <- c(listX, list(lin.d=ssD$X))
		listZ <- c(listZ, list(lin.d=ssD$Z))
	}
	
	if (AA.pop) {
		ssAA <- AxA.effects(dd$Mother_pop, dd$Father_pop, dd$Gen)
		listX <- c(listX, list(pop.aa=ssAA$X))
		listZ <- c(listZ, list(pop.aa=ssAA$Z))
	}
	
	if (AA.line) {
		ssAA <- AxA.effects(dd$Mother_line, dd$Father_line, dd$Gen)
		listX <- c(listX, list(lin.aa=ssAA$X))
		listZ <- c(listZ, list(lin.aa=ssAA$Z))
	}
	
	if (DD.pop) {
		ssDD <- DxD.effects(dd$Mother_pop, dd$Father_pop, dd$Gen)
		listX <- c(listX, list(pop.dd=ssDD$X))
		listZ <- c(listZ, list(pop.dd=ssDD$Z))
	}
	
	if (DD.line) {
		ssDD <- DxD.effects(dd$Mother_line, dd$Father_line, dd$Gen)
		listX <- c(listX, list(lin.dd=ssDD$X))
		listZ <- c(listZ, list(lin.dd=ssDD$Z))
	}
	
	y <- dd[,trait]
	NAy <- is.na(y)
	y <- y[!NAy]
	
	X <- do.call(cbind, listX)
	colnames(X) <- unlist(lapply(names(listX), function(nn) { cc <- colnames(listX[[nn]]); if (is.null(cc)) nn else paste(nn, cc, sep=".")}))
	X <- X[!NAy,,drop=FALSE]
	
	Z <- do.call(cbind, listZ)
	colnames(Z) <- unlist(lapply(names(listZ), function(nn) { cc <- colnames(listZ[[nn]]); if (is.null(cc)) nn else paste(nn, cc, sep=".")}))
	Z <- Z[!NAy,]
	
	startval <- c(
		c(mean(y), rep(0, ncol(X)-1)),
		rep(0, ncol(Z)),
		rep(1, length(listZ)),
		var(y))

	hglm(y=y, X=X, Z=Z, RandC=sapply(listZ, ncol), startval=startval, verbose=FALSE)
}


admb.wrapper <- function(dd, trait="Weight", A.pop=TRUE, A.line=FALSE, D.pop=TRUE, D.line=FALSE, E.pop=TRUE, E.line=FALSE) {
	
	dd <- dd[!is.na(dd[,trait]),]
		
	tt <- dd[,trait]
	common.pop <- factor(c(as.character(dd$Mother_pop), as.character(dd$Father_pop)))
	common.lin <-  factor(c(as.character(dd$Mother_line), as.character(dd$Father_line)))
	newm.pop   <- factor(as.character(dd$Mother_pop), levels=levels(common.pop))
	newf.pop   <- factor(as.character(dd$Father_pop), levels=levels(common.pop))
	newm.lin   <- factor(as.character(dd$Mother_line), levels=levels(common.lin))
	newf.lin   <- factor(as.character(dd$Father_line), levels=levels(common.lin))
	npop <- length(levels(common.pop))
	nlin <- length(levels(common.lin))
	
	admb.data <- list(
		nobs=nrow(dd),
		npop=npop,
		nlin=nlin,
		pop_mother=as.numeric(newm.pop),
		pop_father=as.numeric(newf.pop),
		lin_mother=as.numeric(newm.lin),
		lin_father=as.numeric(newf.lin),
		trait=tt
	)
	
	admb.params <- list(
		mu=mean(tt), 
		mean_d=0,
		mean_e=0,
		sigma_a_pop=sd(tt)/10,
		sigma_a_lin=sd(tt)/10,
		sigma_d_pop=sd(tt)/10,
		sigma_d_lin=sd(tt)/10,
		sigma_e_pop=0.1,
		sigma_e_lin=0.1,
		sigma_R=sd(tt)
	)
	
	admb.re <- list(
		u_a_pop=npop,
		u_a_lin=nlin,
		u_d_pop=npop*(npop-1)/2,
		u_d_lin=nlin*(nlin-1)/2,
		u_e_pop=npop*(npop-1)/2,
		u_e_lin=nlin*-nlin-1).2
	)
	
	my.params <- c("mu","sigma_R")
	my.re <- character(0)
	my.model <- "../src/admb-models/model"
	
	if (A.pop) {
		my.params <- c(my.params, "sigma_a_pop")
		my.re     <- c(my.re,     "u_a_pop")
		my.model <- c(my.model, "aP")
	}
	if (A.line) {
		my.params <- c(my.params, "sigma_a_lin")
		my.re     <- c(my.re,     "u_a_lin")
		my.model <- c(my.model, "aL")		
	}
	if (D.pop) {
		my.params <- c(my.params, "sigma_d_pop", "mean_d")
		my.re     <- c(my.re,     "u_d_pop")
		my.model <- c(my.model, "dP")
	}
	if (D.line) {
		my.params <- c(my.params, "sigma_d_lin", "mean_d")
		my.re     <- c(my.re,     "u_d_lin")
		my.model <- c(my.model, "dL")		
	}
	if (E.pop) {
		my.params <- c(my.params, "sigma_e_pop", "mean_e")
		my.re     <- c(my.re,     "u_e_pop")
		my.model <- c(my.model, "eP")
	}
	if (E.line) {
		my.params <- c(my.params, "sigma_e_lin", "mean_e")
		my.re     <- c(my.re,     "u_e_lin")
		my.model <- c(my.model, "eL")		
	}	
	
	do_admb(paste(my.model, collapse="-"), data=admb.data, param=admb.params[unique(my.params)], re=admb.re[my.re])
}

a.A <- admb.wrapper(dd)

m.A       <- hglm.wrapper(dd, A.pop=TRUE, D.pop=FALSE, AA.pop=FALSE, DD.pop=FALSE)
#~ m.A.D     <- hglm.wrapper(dd, A.pop=TRUE, D.pop=TRUE,  AA.pop=FALSE, DD.pop=FALSE)
#~ m.A.D.AA  <- hglm.wrapper(dd, A.pop=TRUE, D.pop=TRUE,  AA.pop=TRUE,  DD.pop=FALSE)
#~ m.A.D.DD  <- hglm.wrapper(dd, A.pop=TRUE, D.pop=TRUE,  AA.pop=FALSE, DD.pop=TRUE)

#~ # The last one does not converge
#~ m.A.D.AA.DD<-hglm.wrapper(dd, A.pop=TRUE, D.pop=TRUE,  AA.pop=TRUE,  DD.pop=TRUE)


#~ ml.A       <- hglm.wrapper(dd, A.pop=FALSE, D.pop=FALSE, AA.pop=FALSE, DD.pop=FALSE, A.line=TRUE, D.line=FALSE, AA.line=FALSE, DD.line=FALSE)
#~ ml.A.D     <- hglm.wrapper(dd, A.pop=FALSE, D.pop=FALSE, AA.pop=FALSE, DD.pop=FALSE, A.line=TRUE, D.line=TRUE, AA.line=FALSE, DD.line=FALSE)
#~ ml.A.D.AA  <- hglm.wrapper(dd, A.pop=FALSE, D.pop=FALSE, AA.pop=FALSE, DD.pop=FALSE, A.line=TRUE, D.line=TRUE, AA.line=TRUE, DD.line=FALSE)
#~ ml.A.D.DD  <- hglm.wrapper(dd, A.pop=FALSE, D.pop=FALSE, AA.pop=FALSE, DD.pop=FALSE, A.line=TRUE, D.line=TRUE, AA.line=FALSE, DD.line=TRUE)

#~ # Does not converge
#~ ml.A.D.AA.DD<-hglm.wrapper(dd, A.pop=FALSE, D.pop=FALSE, AA.pop=FALSE, DD.pop=FALSE, A.line=TRUE, D.line=TRUE, AA.line=TRUE, DD.line=TRUE)

