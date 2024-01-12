
col <- c(
		ref="darkgray", 
		F1 ="orange", 
		F2 ="green",
		D  ="seagreen3",
		AxA="tomato",
		add="orchid2")
		
data.file <- "../data/data_clean.txt"

col.pops <- setNames(c(RColorBrewer::brewer.pal(3, "Purples"), RColorBrewer::brewer.pal(3, "Greens")), nm=c("SA11","SA16", "SA17", "SA2", "SA3", "SA4"))

weight.name <- "Dry biomass (g)"
silique.name <- "Estimated number of siliques"

fig.width  <- 5.9 # 15 cm in inches
fig.height <- 2.0 #  5 cm in inches
fig.pointsize <- 6
fig.mar    <- c(4,4,2,1)
