source("../src/common.R")
source("../src/model-tools.R")

set.seed(1234)

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
	print(setNames(lapply(fixed.summary.line, function(tt)
		data.frame(
			logLik = round(sapply(tt, logLik), digits=2), 
			df     = sapply(tt, function(x) attributes(logLik(x))$df), 
			AIC    = round(sapply(tt, AIC), digits=2),
			DeltaAIC = round(sapply(tt, AIC) - min(sapply(tt, AIC)), digits=2)
		)), nm=c(weight.name, silique.name)))
sink()

sink("../results/TableS5.txt")
	print(setNames(lapply(fixed.summary.pop, function(tt)
		data.frame(
			logLik = round(sapply(tt, logLik), digits=2), 
			df     = sapply(tt, function(x) attributes(logLik(x))$df), 
			AIC    = round(sapply(tt, AIC), digits=2),
			DeltaAIC = round(sapply(tt, AIC) - min(sapply(tt, AIC)), digits=2)
		)), nm=c(weight.name, silique.name)))
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


pdf("../results/Fig2.pdf", width=fig.width, height=2*fig.height, pointsize=fig.pointsize)
	layout(cbind(1:2, 3:4))
	
	par(cex=1, mar=fig.mar) # layout changes the point size
	
	hist(f.W.d,  breaks=20, xlab="Dominance (D) effect", main=weight.name, xlim=xlim.W)
	abline(v=mean(f.W.d, na.rm=TRUE), col="red", lwd=3)
	hist(f.W.aa, breaks=20, xlab="Epistasis (AxA) effect", main="", xlim=xlim.W)
	abline(v=mean(f.W.aa, na.rm=TRUE), col="red", lwd=3)
	
	hist(f.F.d,  breaks=20, xlab="Dominance (D) effect", main=silique.name, xlim=xlim.F)
	abline(v=mean(f.F.d, na.rm=TRUE), col="red", lwd=3)
	hist(f.F.aa, breaks=20, xlab="Epistasis (AxA) effect", main="", xlim=xlim.F)
	abline(v=mean(f.F.aa, na.rm=TRUE), col="red", lwd=3)
dev.off()


