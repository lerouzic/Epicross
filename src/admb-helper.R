
library(R2admb)

source("../src/genet-model.R")

admb.dir <- "admb"

admb.DATA.common <-
"DATA_SECTION

	init_int nobs
	init_int npop
	init_int nlin
	init_int ninter_pop
	init_int ninter_lin
	
	init_vector pop_mother(1,nobs)
	init_vector pop_father(1,nobs)
	init_vector lin_mother(1,nobs)
	init_vector lin_father(1,nobs)
	init_vector pop_inter(1,nobs)
	init_vector lin_inter(1,nobs)
	
	// For simplicity, model coefficients are provided along with the data
	// It could be recalculated here, but faster/easier/more flexible in R
	init_vector model_a_pop(1,nobs)
	init_vector model_a_lin(1,nobs)
	init_vector model_d_pop(1,nobs)
	init_vector model_d_lin(1,nobs)
	init_vector model_aa_pop(1,nobs)
	init_vector model_aa_lin(1,nobs)
	init_vector model_ad_pop(1,nobs)
	init_vector model_ad_lin(1,nobs)
	init_vector model_dd_pop(1,nobs)
	init_vector model_dd_lin(1,nobs)
	
	init_vector trait(1,nobs)"

admb.PARAMETER.common <- "

PARAMETER_SECTION

	init_number mu(1)
	init_bounded_number sigma_R(0.000001,1000000,1)
	
	vector pred_trait(1,nobs)
	
	objective_function_value f"

admb.PARAMETER.specific.1 <- list(
	a.pop  = "
	init_bounded_number sigma_a_pop(0.000001,1000000,1)",
	a.lin  = "
	init_bounded_number sigma_a_lin(0.000001,1000000,1)",
	d.pop  = "
	init_number mean_d_pop(2)
	init_bounded_number sigma_d_pop(0.000001,1000000,2)",
	d.lin  = "
	init_number mean_d_lin(2)
	init_bounded_number sigma_d_lin(0.000001,1000000,2)",
	e.pop  = "
	init_number mean_e_pop(3)
	init_bounded_number sigma_e_pop(0.000001,1000000,3)",
	e.lin  = "
	init_number mean_e_lin(3)
	init_bounded_number sigma_e_lin(0.000001,1000000,3)"
)

admb.PARAMETER.specific.2 <- list(
	a.pop  = "
	random_effects_vector u_a_pop(1,npop,1)",
	a.lin  = "
	random_effects_vector u_a_lin(1,nlin,1)",
	d.pop  = "
	random_effects_vector u_d_pop(1,ninter_pop,2)",
	d.lin  = "
	random_effects_vector u_d_lin(1,ninter_lin,2)",
	e.pop  = "
	random_effects_vector u_e_pop(1,ninter_pop,3)",
	e.lin  = "
	random_effects_vector u_e_lin(1,ninter_lin,3)"
)

admb.PROCEDURE.common.1 <-"

PROCEDURE_SECTION
	
	f = 0;

	// Prior part for random effects"

admb.PROCEDURE.common.2 <- "

	// Predicted phenotypes
	for (int i = 1; i <= nobs; i++) {
		pred_trait[i] = mu"

admb.PROCEDURE.common.3 <- "
		;
	}
	f -= -nobs*log(sigma_R) - 0.5*norm2((pred_trait-trait)/sigma_R);"

admb.PROCEDURE.specific.1 <- list(
	a.pop  = "
	f -= - 0.5*norm2(u_a_pop);",
	a.lin  = "
	f -= - 0.5*norm2(u_a_lin);",
	d.pop  = "
	f -= - 0.5*norm2(u_d_pop);",
	d.lin  = "
	f -= - 0.5*norm2(u_d_lin);",
	e.pop  = "
	f -= - 0.5*norm2(u_e_pop);",
	e.lin  = "
	f -= - 0.5*norm2(u_e_lin);"
)

admb.PROCEDURE.specific.2 <- list(
	a.pop  = "
		                + model_a_pop(i) *  sigma_a_pop*u_a_pop(pop_mother(i)) 
		                + model_a_pop(i) *  sigma_a_pop*u_a_pop(pop_father(i))",
	a.lin  = "
		                + model_a_lin(i) *  sigma_a_lin*u_a_lin(lin_mother(i))
		                + model_a_lin(i) *  sigma_a_lin*u_a_lin(lin_father(i))",
	d.pop  = "
		                + model_d_pop(i) * (mean_d_pop + sigma_d_pop*u_d_pop(pop_inter(i)))",
	d.lin  = "
		                + model_d_lin(i) * (mean_d_lin + sigma_d_lin*u_d_lin(lin_inter(i)))",
	aa.pop = "
		                + model_aa_pop(i)*  sigma_a_pop*u_a_pop(pop_mother(i))
		                                 *  sigma_a_pop*u_a_pop(pop_father(i)) 
		                                 * (mean_e_pop + sigma_e_pop*u_e_pop(pop_inter(i)))",
	aa.lin = "
		                + model_aa_lin(i)*  sigma_a_lin*u_a_lin(lin_mother(i))
		                                 *  sigma_a_lin*u_a_lin(lin_father(i)) 
		                                 * (mean_e_lin + sigma_e_lin*u_e_lin(lin_inter(i)))",
	ad.pop = "
		                + model_ad_pop(i)*  sigma_a_pop*(u_a_pop(pop_mother(i))+u_a_pop(pop_father(i)))
		                                 * (mean_d_pop + sigma_d_pop*u_d_pop(pop_inter(i)))
		                                 * (mean_e_pop + sigma_e_pop*u_e_pop(pop_inter(i)))",
	ad.lin = "
		                + model_ad_lin(i)*  sigma_a_lin*(u_a_lin(lin_mother(i))+u_a_lin(lin_father(i)))
		                                 * (mean_d_lin + sigma_d_lin*u_d_lin(lin_inter(i)))
		                                 * (mean_e_lin + sigma_e_lin*u_e_lin(lin_inter(i)))",
	dd.pop = "
		                + model_dd_pop(i)* (mean_d_pop + sigma_d_pop*u_d_pop(pop_inter(i)))
		                                 * (mean_d_pop + sigma_d_pop*u_d_pop(pop_inter(i)))
		                                 * (mean_e_pop + sigma_e_pop*u_e_pop(pop_inter(i)))",
	dd.lin = "
		                + model_dd_lin(i)* (mean_d_lin + sigma_d_lin*u_d_lin(lin_inter(i))) 
		                                 * (mean_d_lin + sigma_d_lin*u_d_lin(lin_inter(i)))
		                                 * (mean_e_lin + sigma_e_lin*u_e_lin(lin_inter(i)))"
)

admb.TOP_OF_MAIN_SECTION.common <- "

TOP_OF_MAIN_SECTION

	arrmblsize = 400000000L;
	// gradient_structure::set_ARRAY_MEMBLOCK_SIZE(200000L);
	gradient_structure::set_GRADSTACK_BUFFER_SIZE(100000000L);
	gradient_structure::set_MAX_NVAR_OFFSET(50000);
	gradient_structure::set_CMPDIF_BUFFER_SIZE(50000000L);
"


make.admb <- function(A.pop=TRUE, A.line=FALSE, D.pop=TRUE, D.line=FALSE, E.pop=TRUE, E.line=FALSE) {

	model.name <- paste0("model", 
		if (A.pop)  "-ap" else "", 
		if (A.line) "-al" else "",
		if (D.pop)  "-dp" else "",
		if (D.line) "-dl" else "", 
		if (E.pop)  "-ep" else "",
		if (E.line) "-el" else "")
	model.dir <- file.path(admb.dir, model.name)
	if (!dir.exists(model.dir))
		dir.create(model.dir)
	file.name <- paste0(model.dir, "/", model.name, ".tpl")

	eff <- NULL
	if (A.pop)  eff <- c(eff, "a.pop")
	if (A.line) eff <- c(eff, "a.lin")
	if (D.pop)  eff <- c(eff, "d.pop")
	if (D.line) eff <- c(eff, "d.lin")
	if (E.pop)  eff <- c(eff, "e.pop")
	if (E.line) eff <- c(eff, "e.lin")
	if (A.pop && E.pop)   eff <- c(eff, "aa.pop")
	if (A.line && E.line) eff <- c(eff, "aa.lin")
	if (D.pop && E.pop)   eff <- c(eff, "dd.pop")
	if (D.line && E.line) eff <- c(eff, "dd.lin")
	if (A.pop && D.pop && E.pop)    eff <- c(eff, "ad.pop")
	if (A.line && D.line && E.line) eff <- c(eff, "ad.lin")

	cat(file=file.name, admb.DATA.common)
	cat(file=file.name, admb.PARAMETER.common, append=TRUE)
	cat(file=file.name, do.call(paste, c(admb.PARAMETER.specific.1[eff], list(sep=""))), append=TRUE)
	cat(file=file.name, do.call(paste, c(admb.PARAMETER.specific.2[eff], list(sep=""))), append=TRUE)
	cat(file=file.name, admb.PROCEDURE.common.1, append=TRUE)
	cat(file=file.name, do.call(paste, c(admb.PROCEDURE.specific.1[eff], list(sep=""))), append=TRUE)
	cat(file=file.name, admb.PROCEDURE.common.2, append=TRUE)
	cat(file=file.name, do.call(paste, c(admb.PROCEDURE.specific.2[eff], list(sep=""))), append=TRUE)
	cat(file=file.name, admb.PROCEDURE.common.3, append=TRUE)
	cat(file=file.name, admb.TOP_OF_MAIN_SECTION.common, append=TRUE)

	return(file.path(model.dir, model.name))
}



admb.wrapper <- function(dd, trait="Weight", A.pop=TRUE, A.line=FALSE, D.pop=TRUE, D.line=FALSE, E.pop=TRUE, E.line=FALSE) {
	.ac <- as.character
	
	dd <- dd[!is.na(dd[,trait]),]
	
	tt <- dd[,trait]
	gen <- dd$Gen
	common.pop <- factor(c(.ac(dd$Mother_pop), .ac(dd$Father_pop)))
	common.lin <- factor(c(.ac(dd$Mother_line), .ac(dd$Father_line)))
	newm.pop   <- factor(.ac(dd$Mother_pop), levels=levels(common.pop))
	newf.pop   <- factor(.ac(dd$Father_pop), levels=levels(common.pop))
	newm.lin   <- factor(.ac(dd$Mother_line), levels=levels(common.lin))
	newf.lin   <- factor(.ac(dd$Father_line), levels=levels(common.lin))
	inter.pop  <- factor(ifelse(newm.pop == newf.pop, 
		"self", 
		ifelse(as.numeric(newm.pop) < as.numeric(newf.pop), 
			paste(.ac(newm.pop), .ac(newf.pop), sep="-"), 
			paste(.ac(newf.pop), .ac(newm.pop), sep="-"))))
	inter.lin  <- factor(ifelse(newm.lin == newf.lin, 
		"self",
		ifelse(as.numeric(newm.lin) < as.numeric(newf.lin), 
			paste(.ac(newm.lin), .ac(newf.lin), sep="-"), 
			paste(.ac(newf.lin), .ac(newm.lin), sep="-"))))
	
	# Model coeffcients (makes the code in ADMB easier)
	# The table "MODEL" is in structure.R
	zz <- rep(0, nrow(dd))
	model.a.pop  <- if (A.pop)  ifelse(newm.pop == newf.pop, MODEL["P","a"],  ifelse(gen==1, MODEL["F1","a"],  MODEL["F2","a"]))  else zz
	model.a.lin  <- if (A.line) ifelse(newm.lin == newf.lin, MODEL["P","a"],  ifelse(gen==1, MODEL["F1","a"],  MODEL["F2","a"]))  else zz
	model.d.pop  <- if (D.pop)  ifelse(newm.pop == newf.pop, MODEL["P","d"],  ifelse(gen==1, MODEL["F1","d"],  MODEL["F2","d"]))  else zz
	model.d.lin  <- if (D.line) ifelse(newm.lin == newf.lin, MODEL["P","d"],  ifelse(gen==1, MODEL["F1","d"],  MODEL["F2","d"]))  else zz
	model.aa.pop <- if (E.pop)  ifelse(newm.pop == newf.pop, MODEL["P","aa"], ifelse(gen==1, MODEL["F1","aa"], MODEL["F2","aa"])) else zz
	model.aa.lin <- if (E.line) ifelse(newm.lin == newf.lin, MODEL["P","aa"], ifelse(gen==1, MODEL["F1","aa"], MODEL["F2","aa"])) else zz
	model.ad.pop <- if (E.pop)  ifelse(newm.pop == newf.pop, MODEL["P","ad"], ifelse(gen==1, MODEL["F1","ad"], MODEL["F2","ad"])) else zz
	model.ad.lin <- if (E.line) ifelse(newm.lin == newf.lin, MODEL["P","ad"], ifelse(gen==1, MODEL["F1","ad"], MODEL["F2","ad"])) else zz
	model.dd.pop <- if (E.pop)  ifelse(newm.pop == newf.pop, MODEL["P","dd"], ifelse(gen==1, MODEL["F1","dd"], MODEL["F2","dd"])) else zz
	model.dd.lin <- if (E.line) ifelse(newm.lin == newf.lin, MODEL["P","dd"], ifelse(gen==1, MODEL["F1","dd"], MODEL["F2","dd"])) else zz
	
	npop <- length(levels(common.pop))
	nlin <- length(levels(common.lin))
	ninter.pop <- length(levels(inter.pop))
	ninter.lin <- length(levels(inter.lin))
	
	# Order matters!
	admb.data <- list(
		nobs         = nrow(dd),
		npop         = npop,
		nlin         = nlin,
		ninter_pop   = ninter.pop,
		ninter_lin   = ninter.lin,
		
		pop_mother   = as.numeric(newm.pop),
		pop_father   = as.numeric(newf.pop),
		lin_mother   = as.numeric(newm.lin),
		lin_father   = as.numeric(newf.lin),
		pop_inter    = as.numeric(inter.pop),
		lin_inter    = as.numeric(inter.lin),
		
		model_a_pop  = model.a.pop, 
		model_a_lin  = model.a.lin,
		model_d_pop  = model.d.pop, 
		model_d_lin  = model.d.lin, 
		model_aa_pop = model.aa.pop,
		model_aa_lin = model.aa.lin, 
		model_ad_pop = model.ad.pop, 
		model_ad_lin = model.ad.lin,
		model_dd_pop = model.dd.pop, 
		model_dd_lin = model.dd.lin,
		
		trait        = tt
	)
	
	admb.params <- list(
		mu=mean(tt), 
		mean_d_pop=0,
		mean_d_lin=0,
		mean_e_pop=0,
		mean_e_lin=0,
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
		u_d_pop=ninter.pop,
		u_d_lin=ninter.lin,
		u_e_pop=ninter.pop,
		u_e_lin=ninter.lin
	)
	
	my.params <- c("mu","sigma_R")
	my.re <- character(0)
	
	
	if (A.pop) {
		my.params <- c(my.params, "sigma_a_pop")
		my.re     <- c(my.re,     "u_a_pop")
	}
	if (A.line) {
		my.params <- c(my.params, "sigma_a_lin")
		my.re     <- c(my.re,     "u_a_lin")
	}
	if (D.pop) {
		my.params <- c(my.params, "mean_d_pop", "sigma_d_pop")
		my.re     <- c(my.re,     "u_d_pop")
	}
	if (D.line) {
		my.params <- c(my.params, "mean_d_lin", "sigma_d_lin")
		my.re     <- c(my.re,     "u_d_lin")
	}
	if (E.pop) {
		my.params <- c(my.params, "mean_e_pop", "sigma_e_pop")
		my.re     <- c(my.re,     "u_e_pop")
	}
	if (E.line) {
		my.params <- c(my.params, "mean_e_lin", "sigma_e_lin")
		my.re     <- c(my.re,     "u_e_lin")
	}	
	
	my.model <- make.admb(A.pop, A.line, D.pop, D.line, E.pop, E.line)
	wd <- getwd()
	setwd(dirname(my.model)) # do_admb does not deal well with paths. Let's move into the model directory
	ans <- try(do_admb(basename(my.model), data=admb.data, param=admb.params[my.params], re=admb.re[my.re], verbose=FALSE))
	setwd(wd)
	ans
}


admb.summ <- function(obj, what="estimates") {
	stopifnot("admb" %in% class(obj))
	
	if (what == "estimates") 
		ans <- list(
			mu     = coef(obj)["mu"],
			pop.d  = coef(obj)["mean_d_pop"],
			pop.e  = coef(obj)["mean_e_pop"],
			lin.d  = coef(obj)["mean_d_lin"],
			lin.e  = coef(obj)["mean_e_lin"],
			Vpop.a = coef(obj)["sigma_a_pop"]^2,
			Vlin.a = coef(obj)["sigma_a_lin"]^2,
			Vpop.d = coef(obj)["sigma_d_pop"]^2,
			Vlin.d = coef(obj)["sigma_d_lin"]^2,
			Vpop.e = coef(obj)["sigma_e_pop"]^2,
			Vlin.e = coef(obj)["sigma_e_lin"]^2,
			VR     = coef(obj)["sigma_R"]^2
		)
	if (what == "variances")
		ans <- list(
		)
	ans
}
