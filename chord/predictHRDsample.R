#!/usr/bin/Rscript
library(randomForest)
args <- commandArgs(TRUE)

base_dir <- '/data/common/tools/chord_v0.1/'
setwd(args[1])
hmf_patients_mut_sigs <- read.table(args[3])

rf_model <- readRDS(paste0(base_dir, 'rf_hrd_predict.rds'))

pred <- predict(object = rf_model, newdata = hmf_patients_mut_sigs, type = "prob")
pred <- as.data.frame(pred)

pred$hrd <- pred$BRCA1 + pred$BRCA2

cutoff <- 0.5

pred$predicted_response <- as.integer(pred$hrd >= cutoff)
fileName <- paste(args[2], "_chord_prediction.txt", sep = "")
write.table(pred, fileName, sep = '\t', quote = F)