library(data.table)
library(dplyr)
library(caTools)
library(randomForest)
library(textshape)

#carga de los resultados de GSEA
ddr_gseaResult_neg <- fread("data/CCND1_neg_DDR.xls")
ddr_gseaResult_neg <- ddr_gseaResult_neg[,-c(2,3,4,9,10,11,12)]
ddr_gseaResult_pos <- fread("data/CCND1_pos_DDR.xls")
ddr_gseaResult_pos <- ddr_gseaResult_pos[,-c(2,3,4,9,10,11,12)]

mcl_gseaResult_neg <- fread("data/CCND1_neg_MCL.xls")
mcl_gseaResult_neg <- mcl_gseaResult_neg[,-c(2,3,4,9,10,11,12)]
mcl_gseaResult_pos <- fread("data/CCND1_pos_MCL.xls")
mcl_gseaResult_pos <- mcl_gseaResult_pos[,-c(2,3,4,9,10,11,12)]

cmama_gseaResult_neg <- fread("data/CCND1_neg_CancerDeMama.xls")
cmama_gseaResult_neg <- cmama_gseaResult_neg[,-c(2,3,4,9,10,11,12)]
cmama_gseaResult_pos <- fread("data/CCND1_pos_CancerDeMama.xls")
cmama_gseaResult_pos <- cmama_gseaResult_pos[,-c(2,3,4,9,10,11,12)]

#selección del top 500 de NES
ddr_neg_top500NES <- top_n(ddr_gseaResult_neg, 1000, "NES")[1:1000,]
cmama_neg_top500NES <- top_n(cmama_gseaResult_neg, 1000, "NES")[1:1000,]
mcl_neg_top500NES <- top_n(mcl_gseaResult_neg, 1000, "NES")[1:1000,]
ddr_pos_top500NES <- top_n(ddr_gseaResult_pos, 1000, "NES")[1:1000,]
cmama_pos_top500NES <- top_n(cmama_gseaResult_pos, 1000, "NES")[1:1000,]
mcl_pos_top500NES <- top_n(mcl_gseaResult_pos, 1000, "NES")[1:1000,]

#preparación de ficheros para detectar elementos comunes
fwrite(list(ddr_neg_top500NES$NAME), file="ddr_neg.txt")
fwrite(list(cmama_neg_top500NES$NAME), file="cmama_neg.txt")
fwrite(list(mcl_neg_top500NES$NAME), file="mcl_neg.txt")
fwrite(list(ddr_pos_top500NES$NAME), file="ddr_pos.txt")
fwrite(list(cmama_pos_top500NES$NAME), file="cmama_pos.txt")
fwrite(list(mcl_pos_top500NES$NAME), file="mcl_pos.txt")

#preparacion para random forest
ddr_gseaResult_pos$COMUN <- 0
for (row in 1:nrow(ddr_gseaResult_pos)) {
  if(ddr_gseaResult_pos[row]$NAME %in% cmama_gseaResult_pos$NAME && 
     ddr_gseaResult_pos[row]$NAME %in% mcl_gseaResult_pos$NAME){
    ddr_gseaResult_pos[row]$COMUN <- 1
  }
}
ddr_gseaResult_pos <- column_to_rownames(ddr_gseaResult_pos, "NAME")
names(ddr_gseaResult_pos) <- c("ES", "NES", "p", "FDR", "COMUN")
ddr_gseaResult_pos$COMUN <- as.factor(ddr_gseaResult_pos$COMUN)

ddr_gseaResult_neg$COMUN <- 0
for (row in 1:nrow(ddr_gseaResult_neg)) {
  if(ddr_gseaResult_neg[row]$NAME %in% cmama_gseaResult_neg$NAME && 
     ddr_gseaResult_neg[row]$NAME %in% mcl_gseaResult_neg$NAME){
    ddr_gseaResult_neg[row]$COMUN <- 1
  }
}
ddr_gseaResult_neg <- column_to_rownames(ddr_gseaResult_neg, "NAME")
names(ddr_gseaResult_neg) <- c("ES", "NES", "p", "FDR", "COMUN")
ddr_gseaResult_neg$COMUN <- as.factor(ddr_gseaResult_neg$COMUN)

#aplicación de la técnica
sample_pos = sample.split(ddr_gseaResult_pos$COMUN, SplitRatio = .75)
train_pos = subset(ddr_gseaResult_pos, sample_pos == TRUE)
test_pos  = subset(ddr_gseaResult_pos, sample_pos == FALSE)
rf_pos <- randomForest(COMUN ~ ., data=train_pos)
predicciones_pos <- predict(rf_pos, test_pos)
(mc_pos <- with(test_pos,table(predicciones_pos, COMUN)))
100 * sum(diag(mc_pos)) / sum(mc_pos)

sample_neg = sample.split(ddr_gseaResult_neg$COMUN, SplitRatio = .75)
train_neg = subset(ddr_gseaResult_neg, sample_neg == TRUE)
test_neg  = subset(ddr_gseaResult_neg, sample_neg == FALSE)
rf_neg <- randomForest(COMUN ~ ., data=train_neg)
predicciones_neg <- predict(rf_neg, test_neg)
(mc_neg <- with(test_neg,table(predicciones_neg, COMUN)))
100 * sum(diag(mc_neg)) / sum(mc_neg)


#validación
mtry_pos <- tuneRF(ddr_gseaResult_pos[,-5], ddr_gseaResult_pos$COMUN)
sample_pos = sample.split(ddr_gseaResult_pos$COMUN, SplitRatio = .75)
train_pos = subset(ddr_gseaResult_pos, sample_pos == TRUE)
test_pos  = subset(ddr_gseaResult_pos, sample_pos == FALSE)
rf_pos <- randomForest(COMUN ~ ., data=train_pos, mtry=mtry_pos)
predicciones_pos <- predict(rf_pos, test_pos)
(mc_pos <- with(test_pos,table(predicciones_pos, COMUN)))
100 * sum(diag(mc_pos)) / sum(mc_pos)

mtry_neg <- tuneRF(ddr_gseaResult_neg[,-5], ddr_gseaResult_neg$COMUN)
sample_neg = sample.split(ddr_gseaResult_neg$COMUN, SplitRatio = .75)
train_neg = subset(ddr_gseaResult_neg, sample_neg == TRUE)
test_neg  = subset(ddr_gseaResult_neg, sample_neg == FALSE)
rf_neg <- randomForest(COMUN ~ ., data=train_neg, mtry=mtry_neg)
predicciones_neg <- predict(rf_neg, test_neg)
(mc_neg <- with(test_neg,table(predicciones_neg, COMUN)))
100 * sum(diag(mc_neg)) / sum(mc_neg)

#estratificación del muestreo
mtry_pos <- tuneRF(ddr_gseaResult_pos[,-5], ddr_gseaResult_pos$COMUN)
sample_pos = sample.split(ddr_gseaResult_pos$COMUN, SplitRatio = .75)
train_pos = subset(ddr_gseaResult_pos, sample_pos == TRUE)
test_pos  = subset(ddr_gseaResult_pos, sample_pos == FALSE)
rf_pos <- randomForest(COMUN ~ ., data=train_pos, mtry=mtry_pos, strata = ddr_gseaResult_pos$COMUN, sampsize=c('0'=25, '1'=75))
predicciones_pos <- predict(rf_pos, test_pos)
(mc_pos <- with(test_pos,table(predicciones_pos, COMUN)))
100 * sum(diag(mc_pos)) / sum(mc_pos)

mtry_neg <- tuneRF(ddr_gseaResult_neg[,-5], ddr_gseaResult_neg$COMUN)
sample_neg = sample.split(ddr_gseaResult_neg$COMUN, SplitRatio = .75)
train_neg = subset(ddr_gseaResult_neg, sample_neg == TRUE)
test_neg  = subset(ddr_gseaResult_neg, sample_neg == FALSE)
rf_neg <- randomForest(COMUN ~ ., data=train_neg, mtry=mtry_neg, strata = ddr_gseaResult_neg$COMUN, sampsize=c('0'=25, '1'=75))
predicciones_neg <- predict(rf_neg, test_neg)
(mc_neg <- with(test_neg,table(predicciones_neg, COMUN)))
100 * sum(diag(mc_neg)) / sum(mc_neg)
