#!/usr/bin/Rscript
library(randomForest)
args <- commandArgs(TRUE)

base_dir <- '/data/experiments/hrd_classifier_evaluation/hrdetect_v2/'

hmf_patients_mut_sigs <- read.table(args[1])

rf_model <- readRDS(paste0(base_dir, 'rf_hrd_predict.rds'))

pred <- predict(object = rf_model, newdata = hmf_patients_mut_sigs, type = "prob")
pred <- as.data.frame(pred)

pred$hrd <- pred$BRCA1 + pred$BRCA2

cutoff <- 0.5

pred$predicted_response <- as.integer(pred$hrd >= cutoff)

write.table(pred, 'predictions_sample.txt', sep = '\t', quote = F)