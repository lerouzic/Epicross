
col <- c(
		ref="darkgray", 
		F1 ="orange", 
		F2 ="green",
		D  ="seagreen3",
		AxA="tomato",
		add="orchid2")
		
data.file <- "../data/data_clean.txt"

col.pops <- setNames(rev(viridis::magma(10)), nm=c("SA11","SA16", "SA17", "SA2", "SA3", "SA4"))[1:6]

weight.name <- "Dry Weight (g)"
silique.name <- "Number of siliques"
