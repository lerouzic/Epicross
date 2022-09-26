source("../src/genet-model.R")

library(viridis)

barplot.AkaikeW <- function(delta.AIC, what=c("Weight", "Fitness")) {
	layout(seq_along(what))
	par(mar=c(6, 4, 5, 1))
	for (ww in what) {
		AkW <- exp(-0.5*delta.AIC[[ww]])/sum(exp(-0.5*delta.AIC[[ww]]))
		barplot(AkW, col=magma(length(AkW)), las=1, ylab="Akaike Weight")
		title(ww, line=2.5)
	}
}

plot.fixef <- function(x, jitter=0.05, ...) {
	stopifnot(class(x) == "lm")
	
	RandC <- table(sub("\\..*$", "", sub("Z", "", names(coef(x)[-1])), perl=TRUE))
	
	plot(NULL, xlim=c(0.5, length(RandC)+0.5), ylim=range(x$coef[-1], na.rm=TRUE), xaxt="n", xlab="", ylab="Fixed effect", ...)
	abline(h=0, lty=3)
	axis(1, at=1:length(RandC), labels=names(RandC))
	for (i in seq_along(RandC)) {
		yy <- coef(x)[-1][(c(0,cumsum(RandC))[i]+1):cumsum(RandC)[i]]
		points(rnorm(length(yy), i, sd=jitter), yy)
	}
}

plot.ranef <- function(x, jitter=0.05, ...) {
	stopifnot(class(x) == "hglm")
	
	plot(NULL, xlim=c(0.5, length(x$RandC)+0.5), ylim=range(x$ranef), xaxt="n", xlab="", ylab="Random effect", ...)
	abline(h=0, lty=3)
	axis(1, at=1:length(x$RandC), labels=names(x$RandC))
	for (i in seq_along(x$RandC)) {
		yy <- x$ranef[(c(0,cumsum(x$RandC))[i]+1):cumsum(x$RandC)[i]]
		points(rnorm(length(yy), i, sd=jitter), yy+if(i==1) 0 else x$fixef[i])
	}
}

comp.effects <- function(model1, model2, col=c(a="black", d="red", aa="blue"), xlab="Effect model 1", ylab="Effect model 2", ...) {
	stopifnot(class(model1) %in% c("lm", "hglm"), class(model2) %in% c("lm", "hglm"))
	if (class(model1) == "lm") {
		model1$RandC <- table(sub("\\..*$", "", sub("Z", "", names(coef(model1)[-1])), perl=TRUE))
		model1$ranef <- coef(model1)[-1]
	}
	if (class(model2) == "lm") {
		model2$RandC <- table(sub("\\..*$", "", sub("Z", "", names(coef(model2)[-1])), perl=TRUE))
		model2$ranef <- coef(model2)[-1]
	}
	plot(NULL, xlim=range(model1$ranef, na.rm=TRUE), ylim=range(model2$ranef, na.rm=TRUE), xlab=xlab, ylab=ylab, ...)
	abline(a=0, b=1, lty=3)
	
	for (eff in intersect(names(model1$RandC), names(model2$RandC))) {
		range1 <- (c(0, cumsum(model1$RandC))[which(eff==names(model1$RandC))] +1):cumsum(model1$RandC)[eff]
		range2 <- (c(0, cumsum(model2$RandC))[which(eff==names(model2$RandC))] +1):cumsum(model2$RandC)[eff]
		eff1 <- model1$ranef[range1]
		eff2 <- model2$ranef[range2]
		
		names(eff1) <- sub("^.*\\.", "", names(eff1))
		names(eff2) <- sub("^.*\\.", "", names(eff2))
		
		points(eff1, eff2[names(eff1)], col=col[eff])
	}
}


dd <- read.table("../data/data_clean.txt", stringsAsFactors=FALSE)

# The multilinear analysis does not like "-" in the line names (cannot use the line name as an effect in the model). 
dd$Mother_line <- sub("-", "L", dd$Mother_line)
dd$Father_line <- sub("-", "L", dd$Father_line)


########### Fixed-effect linear model ##############

# Focus on populations

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

fixed.pop.dAIC <- lapply(fixed.pop, function(xx) { aic <- sapply(xx, AIC); aic-min(aic) })

pdf("../results/FnAICWpopLin.pdf", width=15, height=10)
	barplot.AkaikeW(fixed.pop.dAIC)
dev.off()

# Focus on lines

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

fixed.line.dAIC <- lapply(fixed.line, function(xx) { aic <- sapply(xx, AIC); aic-min(aic) })

pdf("../results/FnAICWlineLin.pdf", width=15, height=10)
	barplot.AkaikeW(fixed.line.dAIC)
dev.off()

########### Multilinear model ################

# Focus on populations

fixed.multi.pop <- list(
	Weight = list(
		a            = fit.Npop.Fmu.multi(dd$Mother_pop, dd$Father_pop, dd$Gen, dd$Weight, what=c("a")),
		a.d          = fit.Npop.Fmu.multi(dd$Mother_pop, dd$Father_pop, dd$Gen, dd$Weight, what=c("a", "d")),
		a.aa.e1      = fit.Npop.Fmu.multi(dd$Mother_pop, dd$Father_pop, dd$Gen, dd$Weight, what=c("a", "aa"), e.unique=TRUE),
		a.aa.ee      = fit.Npop.Fmu.multi(dd$Mother_pop, dd$Father_pop, dd$Gen, dd$Weight, what=c("a", "aa"), e.unique=FALSE),
		a.d.aa.e1    = fit.Npop.Fmu.multi(dd$Mother_pop, dd$Father_pop, dd$Gen, dd$Weight, what=c("a", "d", "aa"), e.unique=TRUE),
		a.d.aa.ee    = fit.Npop.Fmu.multi(dd$Mother_pop, dd$Father_pop, dd$Gen, dd$Weight, what=c("a", "d", "aa"), e.unique=FALSE),
		a.d.aa.dd.e1 = fit.Npop.Fmu.multi(dd$Mother_pop, dd$Father_pop, dd$Gen, dd$Weight, what=c("a", "d", "aa", "dd"), e.unique=TRUE),
		a.d.aa.dd.ee = fit.Npop.Fmu.multi(dd$Mother_pop, dd$Father_pop, dd$Gen, dd$Weight, what=c("a", "d", "aa", "dd"), e.unique=FALSE)),
	Fitness = list(
		a            = fit.Npop.Fmu.multi(dd$Mother_pop, dd$Father_pop, dd$Gen, dd$Fitness, what=c("a")),
		a.d          = fit.Npop.Fmu.multi(dd$Mother_pop, dd$Father_pop, dd$Gen, dd$Fitness, what=c("a", "d")),
		a.aa.e1      = fit.Npop.Fmu.multi(dd$Mother_pop, dd$Father_pop, dd$Gen, dd$Fitness, what=c("a", "aa"), e.unique=TRUE),
		a.aa.ee      = fit.Npop.Fmu.multi(dd$Mother_pop, dd$Father_pop, dd$Gen, dd$Fitness, what=c("a", "aa"), e.unique=FALSE),
		a.d.aa.e1    = fit.Npop.Fmu.multi(dd$Mother_pop, dd$Father_pop, dd$Gen, dd$Fitness, what=c("a", "d", "aa"), e.unique=TRUE),
		a.d.aa.ee    = fit.Npop.Fmu.multi(dd$Mother_pop, dd$Father_pop, dd$Gen, dd$Fitness, what=c("a", "d", "aa"), e.unique=FALSE),
		a.d.aa.dd.e1 = fit.Npop.Fmu.multi(dd$Mother_pop, dd$Father_pop, dd$Gen, dd$Fitness, what=c("a", "d", "aa", "dd"), e.unique=TRUE),
		a.d.aa.dd.ee = fit.Npop.Fmu.multi(dd$Mother_pop, dd$Father_pop, dd$Gen, dd$Fitness, what=c("a", "d", "aa", "dd"), e.unique=FALSE)))

fixed.multi.pop.dAIC <- lapply(fixed.multi.pop, function(xx) { aic <- sapply(xx, AIC); aic-min(aic) })

pdf("../results/FnAICWpopMlin.pdf", width=15, height=10)
	barplot.AkaikeW(fixed.multi.pop.dAIC)
dev.off()

# Focus on lines

fixed.multi.line <- list(
	Weight = list(
		a            = fit.Npop.Fmu.multi(dd$Mother_line, dd$Father_line, dd$Gen, dd$Weight, what=c("a")),
		a.d          = fit.Npop.Fmu.multi(dd$Mother_line, dd$Father_line, dd$Gen, dd$Weight, what=c("a", "d")),
		a.aa.e1      = fit.Npop.Fmu.multi(dd$Mother_line, dd$Father_line, dd$Gen, dd$Weight, what=c("a", "aa"), e.unique=TRUE),
		a.aa.ee      = fit.Npop.Fmu.multi(dd$Mother_line, dd$Father_line, dd$Gen, dd$Weight, what=c("a", "aa"), e.unique=FALSE),
		a.d.aa.e1    = fit.Npop.Fmu.multi(dd$Mother_line, dd$Father_line, dd$Gen, dd$Weight, what=c("a", "d", "aa"), e.unique=TRUE),
		a.d.aa.ee    = fit.Npop.Fmu.multi(dd$Mother_line, dd$Father_line, dd$Gen, dd$Weight, what=c("a", "d", "aa"), e.unique=FALSE),
		a.d.aa.dd.e1 = fit.Npop.Fmu.multi(dd$Mother_line, dd$Father_line, dd$Gen, dd$Weight, what=c("a", "d", "aa", "dd"), e.unique=TRUE),
		a.d.aa.dd.ee = fit.Npop.Fmu.multi(dd$Mother_line, dd$Father_line, dd$Gen, dd$Weight, what=c("a", "d", "aa", "dd"), e.unique=FALSE)),
	Fitness = list(
		a            = fit.Npop.Fmu.multi(dd$Mother_line, dd$Father_line, dd$Gen, dd$Fitness, what=c("a")),
		a.d          = fit.Npop.Fmu.multi(dd$Mother_line, dd$Father_line, dd$Gen, dd$Fitness, what=c("a", "d")),
		a.aa.e1      = fit.Npop.Fmu.multi(dd$Mother_line, dd$Father_line, dd$Gen, dd$Fitness, what=c("a", "aa"), e.unique=TRUE),
		a.aa.ee      = fit.Npop.Fmu.multi(dd$Mother_line, dd$Father_line, dd$Gen, dd$Fitness, what=c("a", "aa"), e.unique=FALSE),
		a.d.aa.e1    = fit.Npop.Fmu.multi(dd$Mother_line, dd$Father_line, dd$Gen, dd$Fitness, what=c("a", "d", "aa"), e.unique=TRUE),
		a.d.aa.ee    = fit.Npop.Fmu.multi(dd$Mother_line, dd$Father_line, dd$Gen, dd$Fitness, what=c("a", "d", "aa"), e.unique=FALSE),
		a.d.aa.dd.e1 = fit.Npop.Fmu.multi(dd$Mother_line, dd$Father_line, dd$Gen, dd$Fitness, what=c("a", "d", "aa", "dd"), e.unique=TRUE),
		a.d.aa.dd.ee = fit.Npop.Fmu.multi(dd$Mother_line, dd$Father_line, dd$Gen, dd$Fitness, what=c("a", "d", "aa", "dd"), e.unique=FALSE)))

fixed.multi.line.dAIC <- lapply(fixed.multi.line, function(xx) { aic <- sapply(xx, AIC); aic-min(aic) })

pdf("../results/FnAICWlineMlin.pdf", width=15, height=10)
	barplot.AkaikeW(fixed.multi.line.dAIC)
dev.off()


########### Random-effect linear model ##############

# Focus on populations

rand.pop <- list(
	Weight = list(
		a      = fit.Npop.Fmu.rand(dd$Mother_pop, dd$Father_pop, dd$Gen, dd$Weight, what=c("a")),
		a.d    = fit.Npop.Fmu.rand(dd$Mother_pop, dd$Father_pop, dd$Gen, dd$Weight, what=c("a", "d")),
		a.aa   = fit.Npop.Fmu.rand(dd$Mother_pop, dd$Father_pop, dd$Gen, dd$Weight, what=c("a", "aa")),
		a.d.aa = fit.Npop.Fmu.rand(dd$Mother_pop, dd$Father_pop, dd$Gen, dd$Weight, what=c("a", "d", "aa"))),
	Fitness = list(
		a      = fit.Npop.Fmu.rand(dd$Mother_pop, dd$Father_pop, dd$Gen, dd$Fitness, what=c("a")),
		a.d    = fit.Npop.Fmu.rand(dd$Mother_pop, dd$Father_pop, dd$Gen, dd$Fitness, what=c("a", "d")),
		a.aa   = fit.Npop.Fmu.rand(dd$Mother_pop, dd$Father_pop, dd$Gen, dd$Fitness, what=c("a", "aa")),
		a.d.aa = fit.Npop.Fmu.rand(dd$Mother_pop, dd$Father_pop, dd$Gen, dd$Fitness, what=c("a", "d", "aa"))))

rand.pop.dAIC <- lapply(rand.pop, function(xx) { aic <- sapply(xx, function(xxx) xxx$likelihood$cAIC); aic-min(aic) })
barplot.AkaikeW(rand.pop.dAIC)


# Focus on lines

rand.line <- list(
	Weight = list(
		a      = fit.Npop.Fmu.rand(dd$Mother_line, dd$Father_line, dd$Gen, dd$Weight, what=c("a")),
		a.d    = fit.Npop.Fmu.rand(dd$Mother_line, dd$Father_line, dd$Gen, dd$Weight, what=c("a", "d")),
		a.aa   = fit.Npop.Fmu.rand(dd$Mother_line, dd$Father_line, dd$Gen, dd$Weight, what=c("a", "aa")),
		a.d.aa = fit.Npop.Fmu.rand(dd$Mother_line, dd$Father_line, dd$Gen, dd$Weight, what=c("a", "d", "aa"))),
	Fitness = list(
		a      = fit.Npop.Fmu.rand(dd$Mother_line, dd$Father_line, dd$Gen, dd$Fitness, what=c("a")),
		a.d    = fit.Npop.Fmu.rand(dd$Mother_line, dd$Father_line, dd$Gen, dd$Fitness, what=c("a", "d")),
		a.aa   = fit.Npop.Fmu.rand(dd$Mother_line, dd$Father_line, dd$Gen, dd$Fitness, what=c("a", "aa")),
		a.d.aa = fit.Npop.Fmu.rand(dd$Mother_line, dd$Father_line, dd$Gen, dd$Fitness, what=c("a", "d", "aa"))))

rand.line.dAIC <- lapply(rand.line, function(xx) { aic <- sapply(xx, function(xxx) xxx$likelihood$cAIC); aic-min(aic) })

barplot.AkaikeW(rand.line.dAIC)
