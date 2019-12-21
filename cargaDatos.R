library(data.table) #cargamos la librería para usar data tables
library(tibble) #cargamos la librería para añadir columnas a los data tables
dataDescription_DDR <- fread("data/GPL6947_HumanHT-12_V3_0_R1_11283641_A.bgx", nrows=48803)
data_DDR <- fread("data/GSE25848_series_matrix.txt")
data_DDR <- data_DDR[!is.na(data_DDR$GSM634846)]
data_DDR <- add_column(data_DDR, DESCRIPTION = NA, .after = 1) #añadimos la columna DESCRIPTION con valor NA tal y como indica el manual de uso de GSEA
fwrite(data_DDR, "data/DDR_GSEA.txt", quote=FALSE, sep="\t")

dataDescription_cancerMama <- fread("data/GPL6244-17930.txt", nrows=28869)
data_cancerMama <- fread("data/GSE48989_series_matrix.txt") #cargamos los datos sobre cáncer de mama
data_cancerMama <- data_cancerMama[!is.na(data_cancerMama$GSM1191286)] #eliminamos las filas que tengan valores nulos
data_cancerMama <- add_column(data_cancerMama, DESCRIPTION = NA, .after = 1) #añadimos la columna DESCRIPTION con valor NA tal y como indica el manual de uso de GSEA
fwrite(data_cancerMama, "data/cancerMama_GSEA.txt", quote=FALSE, sep="\t")

dataDescription_MCL <- fread("data/GPL570-55999.txt", nrows=54675)
data_MCL <- fread("data/GSE21452_series_matrix.txt")
data_MCL <- data_MCL[!is.na(data_MCL$GSM536113)]
data_MCL <- add_column(data_MCL, DESCRIPTION = NA, .after = 1) #añadimos la columna DESCRIPTION con valor NA tal y como indica el manual de uso de GSEA
fwrite(data_MCL, "data/MCL_GSEA.txt", quote=FALSE, sep="\t")
