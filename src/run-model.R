source("../src/common.R")
source("../src/model-tools.R")

dd <- read.table(data.file, stringsAsFactors=FALSE)

# The analysis does not like "-" in the line names (cannot use the line name as an effect in the model). 
dd$Mother_line <- sub("-", "L", dd$Mother_line)
dd$Father_line <- sub("-", "L", dd$Father_line)

########### Fixed-effect linear model ##############

fixed.line <- list(
	Weight = list(
		a      = fit.Npop.Fmu(dd$Mother_line, dd$Father_line, dd$Gen, dd$Weight, what=c("a")),
		a.d    = fit.Npop.Fmu(dd$Mother_line, dd$Father_line, dd$Gen, dd$Weight, what=c("a", "d")),
		a.aa   = fit.Npop.Fmu(dd$Mother_line, dd$Father_line, dd$Gen, dd$Weight, what=c("a", "aa")),
		a.d.aa = fit.Npop.Fmu(dd$Mother_line, dd$Father_line, dd$Gen, dd$Weight, what=c("a", "d", "aa"))),
	Fitness = list(
		a      = fit.Npop.Fmu(dd$Mother_line, dd$Father_line, dd$Gen, dd$Fitness, what=c("a")),
		a.d    = fit.Npop.Fmu(dd$Mother_line, dd$Father_line, dd$Gen, dd$Fitness, what=c("a", "d")),
		a.aa   = fit.Npop.Fmu(dd$Mother_line, dd$Father_line, dd$Gen, dd$Fitness, what=c("a", "aa")),
		a.d.aa = fit.Npop.Fmu(dd$Mother_line, dd$Father_line, dd$Gen, dd$Fitness, what=c("a", "d", "aa"))))


fixed.summary.line <- list(
	Weight = list(
		a      = fixed.line$Weight$a, 
		a.d    = fixed.line$Weight$a.d, 
		a.aa   = fixed.line$Weight$a.aa, 
		a.d.aa = fixed.line$Weight$a.d.aa), 
	Fitness = list(
		a      = fixed.line$Fitness$a, 
		a.d    = fixed.line$Fitness$a.d, 
		a.aa   = fixed.line$Fitness$a.aa, 
		a.d.aa = fixed.line$Fitness$a.d.aa)
	)


fixed.pop <- list(
	Weight = list(
		a      = fit.Npop.Fmu(dd$Mother_pop, dd$Father_pop, dd$Gen, dd$Weight, what=c("a")),
		a.d    = fit.Npop.Fmu(dd$Mother_pop, dd$Father_pop, dd$Gen, dd$Weight, what=c("a", "d")),
		a.aa   = fit.Npop.Fmu(dd$Mother_pop, dd$Father_pop, dd$Gen, dd$Weight, what=c("a", "aa")),
		a.d.aa = fit.Npop.Fmu(dd$Mother_pop, dd$Father_pop, dd$Gen, dd$Weight, what=c("a", "d", "aa"))),
	Fitness = list(
		a      = fit.Npop.Fmu(dd$Mother_pop, dd$Father_pop, dd$Gen, dd$Fitness, what=c("a")),
		a.d    = fit.Npop.Fmu(dd$Mother_pop, dd$Father_pop, dd$Gen, dd$Fitness, what=c("a", "d")),
		a.aa   = fit.Npop.Fmu(dd$Mother_pop, dd$Father_pop, dd$Gen, dd$Fitness, what=c("a", "aa")),
		a.d.aa = fit.Npop.Fmu(dd$Mother_pop, dd$Father_pop, dd$Gen, dd$Fitness, what=c("a", "d", "aa"))))


fixed.summary.pop <- list(
	Weight = list(
		a      = fixed.pop$Weight$a, 
		a.d    = fixed.pop$Weight$a.d, 
		a.aa   = fixed.pop$Weight$a.aa, 
		a.d.aa = fixed.pop$Weight$a.d.aa), 
	Fitness = list(
		a      = fixed.pop$Fitness$a, 
		a.d    = fixed.pop$Fitness$a.d, 
		a.aa   = fixed.pop$Fitness$a.aa, 
		a.d.aa = fixed.pop$Fitness$a.d.aa)
	)

################## Writing the AIC tables ################
sink("../results/Table1.txt")
	setNames(lapply(fixed.summary.line, function(tt)
		data.frame(
			logLik = round(sapply(tt, logLik), digits=2), 
			df     = sapply(tt, function(x) attributes(logLik(x))$df), 
			DeltaAIC = round(sapply(tt, AIC) - min(sapply(tt, AIC)), digits=2)
		)), nm=c(weight.name, silique.name))
sink()

sink("../results/TableS4.txt")
	setNames(lapply(fixed.summary.pop, function(tt)
		data.frame(
			logLik = round(sapply(tt, logLik), digits=2), 
			df     = sapply(tt, function(x) attributes(logLik(x))$df), 
			DeltaAIC = round(sapply(tt, AIC) - min(sapply(tt, AIC)), digits=2)
		)), nm=c(weight.name, silique.name))
sink()

################ Distribution of effects #################

f.W.a  <- filter.effects(fixed.line$Weight$a, "a")
f.W.d  <- filter.effects(fixed.line$Weight$a.d, "d")
f.W.aa <- filter.effects(fixed.line$Weight$a.aa, "aa")
f.F.a  <- filter.effects(fixed.line$Fitness$a, "a")
f.F.d  <- filter.effects(fixed.line$Fitness$a.d, "d")
f.F.aa <- filter.effects(fixed.line$Fitness$a.aa, "aa")

xlim.W <- range(c(f.W.d, f.W.aa), na.rm=TRUE)
xlim.F <- range(c(f.F.d, f.F.aa), na.rm=TRUE)


pdf("../results/Fig2.pdf", width=6, height=6)
	layout(cbind(1:2, 3:4))
	
	hist(f.W.d,  breaks=20, xlab="Dominance effect", main=weight.name, xlim=xlim.W)
	abline(v=mean(f.W.d, na.rm=TRUE), col="red", lwd=3)
	hist(f.W.aa, breaks=20, xlab="A x A effect", main="", xlim=xlim.W)
	abline(v=mean(f.W.aa, na.rm=TRUE), col="red", lwd=3)
	
	hist(f.F.d,  breaks=20, xlab="Dominance effect", main=silique.name, xlim=xlim.F)
	abline(v=mean(f.F.d, na.rm=TRUE), col="red", lwd=3)
	hist(f.F.aa, breaks=20, xlab="A x A effect", main="", xlim=xlim.F)
	abline(v=mean(f.F.aa, na.rm=TRUE), col="red", lwd=3)
dev.off()


sink("../results/TableS5.txt")
	setNames(data.frame(
		c(	round(sd  (f.W.a,  na.rm=TRUE), digits=3), 
			round(mean(f.W.d,  na.rm=TRUE), digits=3), 
			round(sd  (f.W.d,  na.rm=TRUE), digits=3), 
			round(mean(f.W.aa, na.rm=TRUE), digits=3), 
			round(sd  (f.W.aa, na.rm=TRUE), digits=3)), 
		c(	round(sd  (f.F.a,  na.rm=TRUE), digits=0), 
			round(mean(f.F.d,  na.rm=TRUE), digits=0), 
			round(sd  (f.F.d,  na.rm=TRUE), digits=0), 
			round(mean(f.F.aa, na.rm=TRUE), digits=0), 
			round(sd  (f.F.aa, na.rm=TRUE), digits=0)), 
		row.names=c("Additive(sd)", "Dominance(mean)", "Dominance(sd)", "Epistasis(mean)" ,"Epistasis(sd)")
		), nm=c(weight.name, silique.name))
sink()
