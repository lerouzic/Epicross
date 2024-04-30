pheno.file    <- "Siliq_count_final.Rform.csv"
position.file <- "Phenotype_data_final_Rform.csv"
cross.file    <- "Croisements.csv"

data.pheno    <- read.csv(pheno.file)
data.position <- read.csv(position.file)
data.cross    <- read.csv(cross.file)

############################ Data fixes and corrections
# in pheno.file
data.pheno$Position[449]  <- "G2_P32_20" # instead of "G2P32_20"
data.pheno$Position[692]  <- "G2_P31_9"  # instead of "G2_31_9"
data.pheno$Position[1061] <- "G2_P32_5"  # instead of "G_P32_5"
data.pheno$Position[1118] <- "G2_P4_22"  # instead of "G2__P4_22"
data.pheno$Position[1121] <- "G2_P4_18"  # instead of "G2__P4_18"

dupl.lines <- c(
    138,             # G1_P1_12  duplicated (lines 12 and 138),   line 138 is mostly NA
     66,             # G1_P25_10 duplicated (lines 66 and 319),   line 66 the weight is very low (35)
    279,             # G1_P23_12 duplicated (lines 279 and 400),  line 179 the weight is very low (37)
    397,             # G1_P23_11 duplicated (lines 397 and 639),  both equivalent, one has to be removed arbitrarily
    645,             # G1_P19_21 duplicated (lines 599 and 645),  line 645 is mostly NA
    419,             # G1_P2_3   duplicated (lines 419 and 651),  line 419 is mostly NA
    925,             # G1_P29_18 duplicated (lines 925 and 949),  both equivalent, one has to be removed arbitrarily
    908,             # G2_P14_10 duplicated (lines 908 and 1012), both equivalent, one has to be removed arbitrarily
    812,             # G2_P14_17 duplicated (lines 812 and 1014), both equivalent, one has to be removed arbitrarily
     43,             # G2_P9_16  duplicated (lines 43 and 1041),  both equivalent, one has to be removed arbitrarily
    198,             # G1_P22_6  duplicated (lines 198 and 1405), line 198 the weight is very low (38)
    394,             # G1_P3_4   duplicated (lines 394 and 1443), both equivalent, one has to be removed arbitrarily
    1009)            # G2_P14_19 duplicated (lines 1009 and 1506),both equivalent, one has to be removed arbitrarily

data.pheno <- data.pheno[-dupl.lines,]     # From this point previous line numbers do not make sense

# To be kept in the data file:
Position <- data.pheno$Position
Mother_line   <- vapply(Position, function(pos) as.character(data.position$Mother[data.position$Position == pos]), FUN.VALUE="A")
Gen           <- vapply(Position, function(pos) as.character(data.position$Gen[data.position$Position == pos]), FUN.VALUE="A")
Cross         <- vapply(Position, function(pos) as.character(data.position$Cross[data.position$Position == pos]),  FUN.VALUE="A")
Father_line   <- vapply(seq_along(Position), function(i) as.character(data.cross$Pere[data.cross$Mere == Mother_line[i] & data.cross$type.cross == Cross[i]]), FUN.VALUE="A")

Weight   <- data.pheno$Weight_tot
Fitness  <- data.pheno$nb_siliq_total_estim

final.data <- data.frame(
    row.names   = Position,
    Mother_line = factor(Mother_line), 
    Mother_pop  = factor(vapply(strsplit(Mother_line, "-"), "[", 1, FUN.VALUE="A")),
    Father_line = factor(Father_line),
    Father_pop  = factor(vapply(strsplit(Father_line, "-"), "[", 1, FUN.VALUE="A")),
    Gen         = factor(Gen),
    Weight      = round(as.numeric(Weight), digits=3),
    Fitness     = round(as.numeric(Fitness), digits=0))

write.table(final.data, file="data_clean.txt", quote=FALSE, sep="\t")
