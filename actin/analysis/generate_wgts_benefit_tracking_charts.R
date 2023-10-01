# Inputs ------------------------------------------------------------------
#change wd to your own path
wd <- paste0(Sys.getenv("HOME"), "/Downloads/")

#add new quartile
Q3_2021 <- c(210701, 210930)
Q4_2021 <- c(211001, 211231)
Q1_2022 <- c(220101, 220331)
Q2_2022 <- c(220401, 220630)
Q3_2022 <- c(220701, 220930)
Q4_2022 <- c(221001, 221231)
Q1_2023 <- c(230101, 230331)
Q2_2023 <- c(230401, 230630)
Q3_2023 <- c(230701, 230930)

#add new quartile
names <- c("2023 Q3", "2023 Q2", "2023 Q1", "2022 Q4", "2022 Q3", "2022 Q2", "2022 Q1", "2021 Q4", "2021 Q3")
quartiles <- list(Q3_2023, Q2_2023, Q1_2023, Q4_2022, Q3_2022, Q2_2022, Q1_2022, Q4_2021, Q3_2021)


# Getting started ---------------------------------------------------------
setwd(wd)
filename_WGS <- "ACTIN WGS+WTS Benefit Tracking - WGS.tsv"
filename_WTS <- "ACTIN WGS+WTS Benefit Tracking - WTS.tsv"
benefit_tracking_WGS <- read.table(file=filename_WGS, sep = '\t', header = TRUE)
benefit_tracking_WTS <- read.table(file=filename_WTS, sep = '\t', header = TRUE)


# LR/CUP ------------------------------------------------------------------
LR_CUP <- data.frame()
i <- 1

for (x in quartiles) {
  LR_CUP[1,i] <- nrow(benefit_tracking_WGS[benefit_tracking_WGS$Category == "LR" & benefit_tracking_WGS$WGS.report.date <= x[2] & benefit_tracking_WGS$WGS.report.date >= x[1], ])
  LR_CUP[2,i] <- nrow(benefit_tracking_WGS[benefit_tracking_WGS$Category == "CUP" & benefit_tracking_WGS$WGS.report.date <= x[2] & benefit_tracking_WGS$WGS.report.date >= x[1], ])
  i <- i+1
}
colnames(LR_CUP) <- names

total_nr_biopsies <- sum(LR_CUP)
total_nr_biopsies_LR <- sum(LR_CUP[1,])
total_nr_biopsies_CUP <- sum(LR_CUP[2,])
total_nr_patients <- subset(benefit_tracking_WGS, benefit_tracking_WGS$WGS.report.date <= quartiles[[1]][2] & benefit_tracking_WGS$WGS.report.date >= quartiles[[length(quartiles)]][1], select = X) 
total_nr_patients_value <- length(unique(total_nr_patients$X))

pdf(file= paste0(wd,"LR_CUP.pdf"), width = 10, height = 7)
par(mar=c(5,5,5,5),mfrow=c(1,1))
barplot(as.matrix(LR_CUP), xlim = c(0,max(LR_CUP[1,])+max(LR_CUP[2,])+5), col=c("blue","red"), las=1, legend = c("LR", "CUP"), horiz = TRUE, args.legend = list(x ='topright', inset = c(-0.05,0.05)), main = paste0("Overall, ", total_nr_biopsies, " biopsies (", total_nr_biopsies_LR, " LR, ", total_nr_biopsies_CUP, " CUP) for ", total_nr_patients_value, " patients have been evaluated in ACTIN "))
grid(nx=NULL,ny=NA,lty=1,col="gray",lwd=1)
barplot(as.matrix(LR_CUP), xlim = c(0,max(LR_CUP[1,])+max(LR_CUP[2,])+5), col=c("blue","red"), las=1, legend = c("LR", "CUP"), horiz = TRUE, args.legend = list(x ='topright', inset = c(-0.05,0.05)), main = paste0("Overall, ", total_nr_biopsies, " biopsies (", total_nr_biopsies_LR, " LR, ", total_nr_biopsies_CUP, " CUP) for ", total_nr_patients_value, " patients have been evaluated in ACTIN "), add=TRUE)
dev.off()


# Sufficient tumor cells --------------------------------------------------
suff_cells <- data.frame()
i <- 1

for (x in quartiles){
  suf <- nrow(benefit_tracking_WGS[benefit_tracking_WGS$WGS.performed.successfully. == "Yes" & benefit_tracking_WGS$WGS.report.date < x[2] & benefit_tracking_WGS$WGS.report.date > x[1], ])
  nonsuf <- nrow(benefit_tracking_WGS[(benefit_tracking_WGS$WGS.performed.successfully. == "No" | benefit_tracking_WGS$WGS.report.date == "")  & benefit_tracking_WGS$WGS.report.date < x[2] & benefit_tracking_WGS$WGS.report.date > x[1], ])
  suff_cells[1,i] <- suf/(suf+nonsuf)*100
  suff_cells[2,i] <- nonsuf/(suf+nonsuf)*100
  
  i <- i+1
}
colnames(suff_cells) <- names

sufficient_nr <- rowMeans(suff_cells, 1)[1]

pdf(file= paste0(wd,"sufficient_tumor_cells.pdf"), width = 10, height = 7)
par(mar=c(5,5,5,5)+.1)
barplot(as.matrix(suff_cells), xlim = c(0,100), col=c("blue","red"), las=1, horiz = TRUE, main = paste0(round(sufficient_nr), "% of biopsies contain sufficient tumor cells for WGTS analysis"), xaxt = "n")
grid(nx=NULL,ny=NA,lty=1,col="gray",lwd=1)
barplot(as.matrix(suff_cells), xlim = c(0,100), col=c("blue","red"), las=1, horiz = TRUE, main = paste0(round(sufficient_nr), "% of biopsies contain sufficient tumor cells for WGTS analysis"), add=TRUE, xaxt = "n")
axis(1,at=c(0,20,40,60,80,100),labels=paste0(c(0,20,40,60,80,100), "%"))
dev.off()


# Sufficient tumor cells LR --------------------------------------------------
suff_cells_LR <- data.frame()
i <- 1

for (x in quartiles){
  suf <- nrow(benefit_tracking_WGS[benefit_tracking_WGS$Category == "LR" & benefit_tracking_WGS$WGS.performed.successfully. == "Yes" & benefit_tracking_WGS$WGS.report.date < x[2] & benefit_tracking_WGS$WGS.report.date > x[1], ])
  nonsuf <- nrow(benefit_tracking_WGS[benefit_tracking_WGS$Category == "LR" & (benefit_tracking_WGS$WGS.performed.successfully. == "No" | benefit_tracking_WGS$WGS.report.date == "")  & benefit_tracking_WGS$WGS.report.date < x[2] & benefit_tracking_WGS$WGS.report.date > x[1], ])
  suff_cells_LR[1,i] <- suf/(suf+nonsuf)*100
  suff_cells_LR[2,i] <- nonsuf/(suf+nonsuf)*100
  
  i <- i+1
}
colnames(suff_cells_LR) <- names

sufficient_nr <- rowMeans(suff_cells_LR, 1)[1]
sufficient_nr_this_quartile <- suff_cells_LR[1, 1]

pdf(file= paste0(wd,"sufficient_tumor_cells_LR.pdf"), width = 10, height = 7)
par(mar=c(5,5,5,5)+.1)
barplot(as.matrix(suff_cells_LR,1), xlim = c(0,100), col=c("blue","red"), las=1, horiz = TRUE, main = "LR", xaxt = "n")
grid(nx=NULL,ny=NA,lty=1,col="gray",lwd=1)
barplot(as.matrix(suff_cells_LR,1), xlim = c(0,100), col=c("blue","red"), las=1, horiz = TRUE, main = "LR", add=TRUE, xaxt = "n")
axis(1,at=c(0,20,40,60,80,100),labels=paste0(c(0,20,40,60,80,100), "%"))
dev.off()


# Sufficient tumor cells CUP --------------------------------------------------
suff_cells_CUP <- data.frame()
i <- 1

for (x in quartiles){
  suf <- nrow(benefit_tracking_WGS[benefit_tracking_WGS$Category == "CUP" & benefit_tracking_WGS$WGS.performed.successfully. == "Yes" & benefit_tracking_WGS$WGS.report.date < x[2] & benefit_tracking_WGS$WGS.report.date > x[1], ])
  nonsuf <- nrow(benefit_tracking_WGS[benefit_tracking_WGS$Category == "CUP" & (benefit_tracking_WGS$WGS.performed.successfully. == "No" | benefit_tracking_WGS$WGS.report.date == "")  & benefit_tracking_WGS$WGS.report.date < x[2] & benefit_tracking_WGS$WGS.report.date > x[1], ])
  suff_cells_CUP[1,i] <- suf/(suf+nonsuf)*100
  suff_cells_CUP[2,i] <- nonsuf/(suf+nonsuf)*100
  
  i <- i+1
}
colnames(suff_cells_CUP) <- names

sufficient_nr <- rowMeans(suff_cells_CUP, 1)[1]

pdf(file= paste0(wd,"sufficient_tumor_cells_CUP.pdf"), width = 10, height = 7)
par(mar=c(5,5,5,5)+.1)
barplot(as.matrix(suff_cells_CUP,1), xlim = c(0,100), col=c("blue","red"), las=1, horiz = TRUE, main = "CUP", xaxt = "n")
grid(nx=NULL,ny=NA,lty=1,col="gray",lwd=1)
barplot(as.matrix(suff_cells_CUP,1), xlim = c(0,100), col=c("blue","red"), las=1, horiz = TRUE, main = "CUP", add=TRUE, xaxt = "n")
axis(1,at=c(0,20,40,60,80,100),labels=paste0(c(0,20,40,60,80,100), "%"))
dev.off()


# TAT in calender days ---------------------------------------------------------------------
tat_calender <- data.frame()
i <- 1

for (x in quartiles) {
  tat <- subset(benefit_tracking_WGS, benefit_tracking_WGS$WGS.report.date < x[2] & benefit_tracking_WGS$WGS.report.date > x[1], select=TAT...calendar.days.)
  tat_calender[1,i] <- mean(tat$TAT...calendar.days.)
  i <- i+1
}
colnames(tat_calender) <- names

average_tat <- mean(benefit_tracking_WGS$TAT...calendar.days.)

pdf(file= paste0(wd,"TAT_calender_days.pdf"), width = 10, height = 7)
barplot(as.matrix(rev(tat_calender)), ylab = "TAT calender days", ylim = c(0,max(tat_calender)+1), col = "blue", main = paste0("Average TAT (turnaround-time) in calendar days for WGS is ", average_tat, " days"))
grid(nx=NA,ny=NULL,lty=1,col="gray",lwd=1)
barplot(as.matrix(rev(tat_calender)), ylab = "TAT calender days", ylim = c(0,max(tat_calender)+1), col = "blue", main = paste0("Average TAT (turnaround-time) in calendar days for WGS is ", average_tat, " days"), add=TRUE)
dev.off()

tat_calender_this_quartile <- tat_calender[1, 1]


# WGS allowed therapy ----------------------------------------------
allowed_therapy <- nrow(benefit_tracking_WGS[benefit_tracking_WGS$WGS.allowed.therapy. == "Yes" & benefit_tracking_WGS$Category == "LR" , ])
non_allowed_therapy <- nrow(benefit_tracking_WGS[benefit_tracking_WGS$WGS.allowed.therapy. == "No" & benefit_tracking_WGS$Category == "LR", ])
slices <- c(allowed_therapy, non_allowed_therapy)

lbls <- c("Yes:", "No:")
pct <- round(slices/sum(slices)*100,1)
lbls <- paste(lbls, pct)
lbls <- paste0(lbls,"%")

pdf(file= paste0(wd,"WGS_allowed_therapy.pdf"), width = 10, height = 7)
pie(slices, labels = lbls, col=c("red","blue"),main="Count of WGS allowed therapy?")
dev.off()


# Patient got treated based on WGS biomarker ----------------------------------------------
yes_treated_based_on_biomarker <- nrow(benefit_tracking_WGS[benefit_tracking_WGS$Patient.got.treated.based.on.WGS.biomarker. == "Yes" & benefit_tracking_WGS$WGS.allowed.therapy. == "Yes" & benefit_tracking_WGS$Category == "LR", ])
no_treated_based_on_biomarker <- nrow(benefit_tracking_WGS[benefit_tracking_WGS$Patient.got.treated.based.on.WGS.biomarker. == "No" & benefit_tracking_WGS$WGS.allowed.therapy. == "Yes" & benefit_tracking_WGS$Category == "LR", ])
unknown_treated_based_on_biomarker <- nrow(benefit_tracking_WGS[benefit_tracking_WGS$Patient.got.treated.based.on.WGS.biomarker. == "Unknown" & benefit_tracking_WGS$WGS.allowed.therapy. == "Yes" & benefit_tracking_WGS$Category == "LR", ])
question_mark_treated_based_on_biomarker <- nrow(benefit_tracking_WGS[benefit_tracking_WGS$Patient.got.treated.based.on.WGS.biomarker. == "?" & benefit_tracking_WGS$WGS.allowed.therapy. == "Yes" & benefit_tracking_WGS$Category == "LR", ])
slices <- c(yes_treated_based_on_biomarker, no_treated_based_on_biomarker, unknown_treated_based_on_biomarker, question_mark_treated_based_on_biomarker)

lbls <- c("Yes:", "No:", "Unknown:", "?:")
pct <- round(slices/sum(slices)*100,1)
lbls <- paste(lbls, pct)
lbls <- paste0(lbls,"%")

pdf(file= paste0(wd,"treated_based_on_WGS.pdf"), width = 10, height = 7)
pie(slices,labels = lbls, col=c("lightgreen","blue","orange","red"),main="Patient got treated based on WGS biomarker?")
dev.off()


# WGS impact category (MIGHT NEED ADJUSTMENTS IF NEW CATEGORIES ARE ADDED)  ----------------------------------------------
impact_category <- data.frame()

no_slots <- nrow(benefit_tracking_WGS[benefit_tracking_WGS$WGS.impact.category == "No slots available for WGS-informed treatment" & benefit_tracking_WGS$WGS.allowed.therapy. == "Yes", ])
non_WGS_pref <- nrow(benefit_tracking_WGS[(benefit_tracking_WGS$WGS.impact.category == "Non-WGS-informed treatment considered preferred" | benefit_tracking_WGS$WGS.impact.category == "Non-WGS-informed treatment considered preferred (screen failure)") & benefit_tracking_WGS$WGS.allowed.therapy. == "Yes", ])
not_meet_req <- nrow(benefit_tracking_WGS[benefit_tracking_WGS$WGS.impact.category == "Patient did not meet non-WGS study requirements" & benefit_tracking_WGS$WGS.allowed.therapy. == "Yes", ])
not_fit <- nrow(benefit_tracking_WGS[benefit_tracking_WGS$WGS.impact.category == "Patient not fit enough for treatment" & benefit_tracking_WGS$WGS.allowed.therapy. == "Yes", ])
not_willing <- nrow(benefit_tracking_WGS[(benefit_tracking_WGS$WGS.impact.category == "Patient not willing to take treatment" | benefit_tracking_WGS$WGS.impact.category == "Patient decided not to want experimental treatment" | benefit_tracking_WGS$WGS.impact.category == "Patient decided not to want to start experimental treatment / not fit enough" | benefit_tracking_WGS$WGS.impact.category == "Patient did not want to receive treatment based off unravelled tumor type" ) & benefit_tracking_WGS$WGS.allowed.therapy. == "Yes", ])
soc_new_diag <- nrow(benefit_tracking_WGS[(benefit_tracking_WGS$WGS.impact.category == "Treated using SOC for uncovered tumor type" | benefit_tracking_WGS$WGS.impact.category == "Treated using SOC for new diagnosis based off WGS" | benefit_tracking_WGS$WGS.impact.category == "Treased using SOC for uncovered tumor type" ) & benefit_tracking_WGS$WGS.allowed.therapy. == "Yes", ])
not_effective <- nrow(benefit_tracking_WGS[benefit_tracking_WGS$WGS.impact.category == "WGS-informed treatment not expected to be effective" & benefit_tracking_WGS$WGS.allowed.therapy. == "Yes", ])

impact_category[1,1] <- no_slots
impact_category[1,2] <- non_WGS_pref
impact_category[1,3] <- not_meet_req
impact_category[1,4] <- not_fit
impact_category[1,5] <- not_willing
impact_category[1,6] <- soc_new_diag
impact_category[1,7] <- not_effective

colnames(impact_category) <- c("No slots available for WGS-informed treatment", "Non-WGS-informed treatment considered preferred", "Patient did not meet non-WGS study requirements", "Patient not fit enough for treatment", "Patient not willing to take treatment", "Treated using SOC for uncovered tumor type", "WGS-informed treatment not expected to be effective")

pdf(file= paste0(wd,"WGS_impact_category.pdf"), width = 20, height = 7)
par(mar=c(5,25,5,5))
barplot(as.matrix(rev(impact_category)), xlim = c(0,35), col = "blue", xlab = "Count of WGS impact category", las=1, horiz = TRUE, main = "WGS impact category if patient not treated based on WGS biomarker")
grid(nx=NULL,ny=NA,lty=1,col="gray",lwd=1)
barplot(as.matrix(rev(impact_category)), xlim = c(0,35), col = "blue", xlab = "Count of WGS impact category", las=1, horiz = TRUE, main = "WGS impact category if patient not treated based on WGS biomarker", add=TRUE)
dev.off()

lr_with_biomarker <- nrow(benefit_tracking_WGS[benefit_tracking_WGS$WGS.allowed.therapy. == "Yes" & benefit_tracking_WGS$Category == "LR", ])
lr_with_biomarker_not_treat <- nrow(benefit_tracking_WGS[benefit_tracking_WGS$WGS.allowed.therapy. == "Yes" & benefit_tracking_WGS$Category == "LR" & benefit_tracking_WGS$Patient.got.treated.based.on.WGS.biomarker. == "No", ])
lr_with_biomarker_not_treated_perc <- (lr_with_biomarker_not_treat / lr_with_biomarker) * 100
lr_with_biomarker_not_treated_not_fit_perc <- (nrow(benefit_tracking_WGS[benefit_tracking_WGS$WGS.allowed.therapy. == "Yes" & benefit_tracking_WGS$Category == "LR" & benefit_tracking_WGS$Patient.got.treated.based.on.WGS.biomarker. == "No" & benefit_tracking_WGS$WGS.impact.category == "Patient not fit enough for treatment", ]) / lr_with_biomarker_not_treat) * 100



# WGS for CUP patients ----------------------------------------------
WGS_did_not_impacted_diagnosis <- nrow(benefit_tracking_WGS[benefit_tracking_WGS$WGS.led.to..change.of..diagnosis. == "No" & benefit_tracking_WGS$Category == "CUP", ])
biopsy_contained_insuff <- nrow(benefit_tracking_WGS[benefit_tracking_WGS$WGS.performed.successfully. == "No" & benefit_tracking_WGS$Category == "CUP", ])
confirmation <- nrow(benefit_tracking_WGS[benefit_tracking_WGS$WGS.led.to.confirmation.of.uncertain.diagnosis. == "Yes" & benefit_tracking_WGS$Category == "CUP", ])
actual_diagnosis <- nrow(benefit_tracking_WGS[benefit_tracking_WGS$WGS.led.to..change.of..diagnosis. == "Yes" & benefit_tracking_WGS$Category == "CUP", ])
confirmation_or_change_diagnosis <- nrow(benefit_tracking_WGS[(benefit_tracking_WGS$WGS.led.to.confirmation.of.uncertain.diagnosis. == "Yes" | benefit_tracking_WGS$WGS.led.to..change.of..diagnosis. == "Yes") & benefit_tracking_WGS$Category == "CUP", ])
slices <- c(WGS_did_not_impacted_diagnosis, biopsy_contained_insuff, confirmation_or_change_diagnosis)

nr_confirmation_or_diagnosis <- round(confirmation_or_change_diagnosis / sum(slices)*100)
nr_confirmation <- round(confirmation / (confirmation + actual_diagnosis)*100)
nr_actual_diagnosis <- round(actual_diagnosis / (confirmation + actual_diagnosis)*100)

lbls <- c("WGS did not impact diagnosis:", "Biopsy contained insufficient tumor cells:", "WGS led to confirmation or change of diagnosis:")
pct <- round(slices/sum(slices)*100,1)
lbls <- paste(lbls, pct)
lbls <- paste0(lbls,"%")

pdf(file= paste0(wd,"WGS_CUP.pdf"), width = 10, height = 7)
par(mar=c(5,7,5,5))
pie(slices,labels = lbls, col=c("orange","blue","red"),main = paste0("For ", nr_confirmation_or_diagnosis, "% CUP patients WGS led to confirmation (", nr_confirmation, "%) or actual diagnosis (", nr_actual_diagnosis, "%)"))
dev.off()


# WTS CUP ----------------------------------------------
WTS_CUP_improved_prediction <- nrow(benefit_tracking_WTS[benefit_tracking_WTS$WTS.improved.CUPPA.prediction.for.determining.primary.tumor.location. == "Yes" & benefit_tracking_WTS$Category == "CUP", ])
nr_of_cups <- nrow(benefit_tracking_WTS[benefit_tracking_WTS$RNA.prep.quality.sufficient. == "Yes" & benefit_tracking_WTS$Category == "CUP", ])

WTS_CUP <- data.frame()
WTS_CUP[1,1] <- WTS_CUP_improved_prediction / nr_of_cups * 100
WTS_CUP[2,1] <- 100 - WTS_CUP_improved_prediction / nr_of_cups * 100
colnames(WTS_CUP) <- ""

pdf(file=paste0(wd,"WTS_CUP.pdf"), width = 10, height = 3)
barplot(as.matrix(WTS_CUP), col=c("blue","red"), las=1, horiz = TRUE, legend = c("Yes", "No"), args.legend = list(x ='right', inset = c(1,1)), xaxt = "n", main = "CUPs where WTS leads to high confidence CUPPA prediction")
grid(nx=NULL,ny=NA,lty=1,col="gray",lwd=1)
barplot(as.matrix(WTS_CUP), col=c("blue","red"), las=1, horiz = TRUE, legend = c("Yes", "No"), args.legend = list(x ='right', inset = c(1,1)), main = "CUPs where WTS leads to high confidence CUPPA prediction", add=TRUE, xaxt = "n")
axis(1,at=c(0,20,40,60,80,100),labels=paste0(c(0,20,40,60,80,100), "%"))
dev.off()


# WTS LR ----------------------------------------------
WTS_LR_potential_actionable <- nrow(benefit_tracking_WTS[benefit_tracking_WTS$WTS.detected.novel.actionable.events..not.detected.or.detectable.by.WGS. == "Yes" & benefit_tracking_WTS$Category == "LR", ])
nr_of_lr <- nrow(benefit_tracking_WTS[benefit_tracking_WTS$RNA.prep.quality.sufficient. == "Yes" & benefit_tracking_WTS$Category == "LR", ])

WTS_LR <- data.frame()
WTS_LR[1,1] <- WTS_LR_potential_actionable / nr_of_lr * 100
WTS_LR[2,1] <- 100 - WTS_LR_potential_actionable / nr_of_lr * 100
colnames(WTS_LR) <- ""

pdf(file=paste0(wd,"WTS_LR.pdf"), width = 10, height = 3)
barplot(as.matrix(WTS_LR), col=c("blue","red"), las=1, horiz = TRUE, legend = c("Yes", "No"), args.legend = list(x ='right', inset = c(1,1)), main = "LR patients with potential additional biomarker from WTS", xaxt = "n")
grid(nx=NULL,ny=NA,lty=1,col="gray",lwd=1)
barplot(as.matrix(WTS_LR), col=c("blue","red"), las=1, horiz = TRUE, legend = c("Yes", "No"), args.legend = list(x ='right', inset = c(1,1)), main = "LR patients with potential additional biomarker from WTS", add=TRUE, xaxt = "n")
axis(1,at=c(0,20,40,60,80,100),labels=paste0(c(0,20,40,60,80,100), "%"))
dev.off()

WTS_succeeded <- nrow(benefit_tracking_WTS[benefit_tracking_WTS$RNA.prep.quality.sufficient. == "Yes", ])
WTS_not_succeeded <- nrow(benefit_tracking_WTS[benefit_tracking_WTS$RNA.prep.quality.sufficient. == "No", ])
WTS_succeeded_perc <- (WTS_succeeded / (WTS_not_succeeded + WTS_succeeded))*100
WTS_hypoth_benefit <- nrow(benefit_tracking_WTS[benefit_tracking_WTS$RNA.prep.quality.sufficient. == "Yes" & (benefit_tracking_WTS$WTS.detected.novel.actionable.events..not.detected.or.detectable.by.WGS. == "Yes" | benefit_tracking_WTS$WTS.improved.CUPPA.prediction.for.determining.primary.tumor.location. == "Yes") , ]) / WTS_succeeded * 100


