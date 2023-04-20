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

################## Writing the AIC table ################
sink("../results/Table1.txt")
	lapply(fixed.summary.line, function(tt)
		data.frame(
			logLik = sapply(tt, logLik), 
			df     = sapply(tt, function(x) attributes(logLik(x))$df), 
			DeltaAIC = sapply(tt, AIC) - min(sapply(tt, AIC))
		)) 
sink()

################ Distribution of effects #################

f.W.d  <- filter.effects(fixed.line$Weight$a.d, "d")
f.W.aa <- filter.effects(fixed.line$Weight$a.aa, "aa")
f.F.d  <- filter.effects(fixed.line$Fitness$a.d, "d")
f.F.aa <- filter.effects(fixed.line$Fitness$a.aa, "aa")

xlim.W <- range(c(f.W.d, f.W.aa), na.rm=TRUE)
xlim.F <- range(c(f.F.d, f.F.aa), na.rm=TRUE)


pdf("../results/Fig2.pdf", width=6, height=6)
	layout(cbind(1:2, 3:4))
	
	hist(f.W.d,  breaks=20, xlab="Dominance effect", main="Weight", xlim=xlim.W)
	abline(v=mean(f.W.d, na.rm=TRUE), col="red", lwd=3)
	hist(f.W.aa, breaks=20, xlab="A x A effect", main="", xlim=xlim.W)
	abline(v=mean(f.W.aa, na.rm=TRUE), col="red", lwd=3)
	
	hist(f.F.d,  breaks=20, xlab="Dominance effect", main="Siliques", xlim=xlim.F)
	abline(v=mean(f.F.d, na.rm=TRUE), col="red", lwd=3)
	hist(f.F.aa, breaks=20, xlab="A x A effect", main="", xlim=xlim.F)
	abline(v=mean(f.F.aa, na.rm=TRUE), col="red", lwd=3)
dev.off()
