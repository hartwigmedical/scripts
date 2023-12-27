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
Q4_2023 <- c(231001, 231231)

#add new quartile
names <- c("2023 Q4", "2023 Q3", "2023 Q2", "2023 Q1", "2022 Q4", "2022 Q3", "2022 Q2", "2022 Q1", "2021 Q4", "2021 Q3")
quartiles <- list(Q4_2023, Q3_2023, Q2_2023, Q1_2023, Q4_2022, Q3_2022, Q2_2022, Q1_2022, Q4_2021, Q3_2021)


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
  LR_CUP[1,i] <- nrow(benefit_tracking_WGS[benefit_tracking_WGS$Category == "LR" & benefit_tracking_WGS$WGS.report.date <= x[2] & benefit_tracking_WGS$WGS.report.date >= x[1] & !is.nan(benefit_tracking_WGS$WGS.report.date), ])
  LR_CUP[2,i] <- nrow(benefit_tracking_WGS[benefit_tracking_WGS$Category == "CUP" & benefit_tracking_WGS$WGS.report.date <= x[2] & benefit_tracking_WGS$WGS.report.date >= x[1] & !is.nan(benefit_tracking_WGS$WGS.report.date), ])
  i <- i+1
}
colnames(LR_CUP) <- names

total_nr_biopsies <- sum(LR_CUP)
total_nr_biopsies_LR <- sum(LR_CUP[1,])
total_nr_biopsies_CUP <- sum(LR_CUP[2,])
total_nr_patients <- subset(benefit_tracking_WGS, benefit_tracking_WGS$WGS.report.date <= quartiles[[1]][2] & benefit_tracking_WGS$WGS.report.date >= quartiles[[length(quartiles)]][1] & !is.nan(benefit_tracking_WGS$WGS.report.date), select = X)
total_nr_patients_value <- length(unique(total_nr_patients$X))

pdf(file= paste0(wd,"LR_CUP.pdf"), width = 10, height = 7)
par(mar=c(5,5,5,5),mfrow=c(1,1))
barplot(as.matrix(LR_CUP), xlim = c(0,max(LR_CUP[1,])+max(LR_CUP[2,])+5), col=c("blue","red"), las=1, legend = c("LR", "CUP"), horiz = TRUE, args.legend = list(x ='topright', inset = c(-0.05,0.05)))
grid(nx=NULL,ny=NA,lty=1,col="gray",lwd=1)
barplot(as.matrix(LR_CUP), xlim = c(0,max(LR_CUP[1,])+max(LR_CUP[2,])+5), col=c("blue","red"), las=1, legend = c("LR", "CUP"), horiz = TRUE, args.legend = list(x ='topright', inset = c(-0.05,0.05)), add=TRUE)
invisible(dev.off())

noquote(paste0("Overall, ", total_nr_biopsies, " biopsies (", total_nr_biopsies_LR, " LR, ", total_nr_biopsies_CUP, " CUP) for ", total_nr_patients_value, " patients have been evaluated in ACTIN "))
noquote(paste0("In ", names[1], " ", LR_CUP[1,1] + LR_CUP[2,1], " WGTS analyses were done (", LR_CUP[1,1], " LR and ", LR_CUP[2,1], " CUP)"))

# Sufficient tumor cells --------------------------------------------------
suff_cells <- data.frame()
i <- 1

for (x in quartiles){
  suf <- nrow(benefit_tracking_WGS[benefit_tracking_WGS$WGS.performed.successfully. == "Yes" & benefit_tracking_WGS$WGS.report.date <= x[2] & benefit_tracking_WGS$WGS.report.date >= x[1] & !is.nan(benefit_tracking_WGS$WGS.report.date), ])
  nonsuf <- nrow(benefit_tracking_WGS[benefit_tracking_WGS$WGS.performed.successfully. == "No" & benefit_tracking_WGS$WGS.report.date <= x[2] & benefit_tracking_WGS$WGS.report.date >= x[1] & !is.nan(benefit_tracking_WGS$WGS.report.date), ])
  suff_cells[1,i] <- suf/(suf+nonsuf)*100
  suff_cells[2,i] <- nonsuf/(suf+nonsuf)*100

  i <- i+1
}
colnames(suff_cells) <- names

sufficient_nr <- rowMeans(suff_cells, 1)[1]

pdf(file= paste0(wd,"sufficient_tumor_cells.pdf"), width = 10, height = 7)
par(mar=c(5,5,5,5)+.1)
barplot(as.matrix(suff_cells), xlim = c(0,100), col=c("blue","red"), las=1, horiz = TRUE, xaxt = "n")
grid(nx=NULL,ny=NA,lty=1,col="gray",lwd=1)
barplot(as.matrix(suff_cells), xlim = c(0,100), col=c("blue","red"), las=1, horiz = TRUE, add=TRUE, xaxt = "n")
axis(1,at=c(0,20,40,60,80,100),labels=paste0(c(0,20,40,60,80,100), "%"))
invisible(dev.off())

noquote(paste0(round(sufficient_nr), "% of biopsies contain sufficient tumor cells for WGTS analysis"))
noquote(paste0("In ", names[1], ", ", round(suff_cells[1,1]), "% of biopsies contained sufficient tumor cells"))

# Sufficient tumor cells LR --------------------------------------------------
suff_cells_LR <- data.frame()
i <- 1

for (x in quartiles){
  suf <- nrow(benefit_tracking_WGS[benefit_tracking_WGS$Category == "LR" & benefit_tracking_WGS$WGS.performed.successfully. == "Yes" & benefit_tracking_WGS$WGS.report.date <= x[2] & benefit_tracking_WGS$WGS.report.date >= x[1] & !is.nan(benefit_tracking_WGS$WGS.report.date), ])
  nonsuf <- nrow(benefit_tracking_WGS[benefit_tracking_WGS$Category == "LR" & benefit_tracking_WGS$WGS.performed.successfully. == "No"  & benefit_tracking_WGS$WGS.report.date <= x[2] & benefit_tracking_WGS$WGS.report.date >= x[1] & !is.nan(benefit_tracking_WGS$WGS.report.date), ])
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
invisible(dev.off())

# Sufficient tumor cells CUP --------------------------------------------------
suff_cells_CUP <- data.frame()
i <- 1

for (x in quartiles){
  suf <- nrow(benefit_tracking_WGS[benefit_tracking_WGS$Category == "CUP" & benefit_tracking_WGS$WGS.performed.successfully. == "Yes" & benefit_tracking_WGS$WGS.report.date <= x[2] & benefit_tracking_WGS$WGS.report.date >= x[1] & !is.nan(benefit_tracking_WGS$WGS.report.date), ])
  nonsuf <- nrow(benefit_tracking_WGS[benefit_tracking_WGS$Category == "CUP" & benefit_tracking_WGS$WGS.performed.successfully. == "No" & benefit_tracking_WGS$WGS.report.date <= x[2] & benefit_tracking_WGS$WGS.report.date >= x[1] & !is.nan(benefit_tracking_WGS$WGS.report.date), ])
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
invisible(dev.off())

noquote(paste0(round(rowMeans(suff_cells_LR)[1]), "% of LR biopsies contain sufficient tumor cells while for CUPs ", round(rowMeans(suff_cells_CUP)[1]), "% of biopsies are sufficient"))

# TAT in calender days ---------------------------------------------------------------------
tat_calender <- data.frame()
i <- 1

for (x in quartiles) {
  tat <- subset(benefit_tracking_WGS, benefit_tracking_WGS$WGS.report.date <= x[2] & benefit_tracking_WGS$WGS.report.date >= x[1] & !is.nan(benefit_tracking_WGS$WGS.report.date), select=TAT...calendar.days.)
  tat_calender[1,i] <- mean(tat$TAT...calendar.days.)
  i <- i+1
}
colnames(tat_calender) <- names

average_tat <- mean(benefit_tracking_WGS$TAT...calendar.days.)

pdf(file= paste0(wd,"TAT_calender_days.pdf"), width = 10, height = 7)
barplot(as.matrix(rev(tat_calender)), ylab = "TAT calender days", ylim = c(0,max(tat_calender)+1), col = "blue")
grid(nx=NA,ny=NULL,lty=1,col="gray",lwd=1)
barplot(as.matrix(rev(tat_calender)), ylab = "TAT calender days", ylim = c(0,max(tat_calender)+1), col = "blue", add=TRUE)
invisible(dev.off())

tat_calender_this_quartile <- tat_calender[1, 1]

noquote(paste0("Average TAT (turnaround-time) in calendar days for WGS is ", round(average_tat,1), " days"))
noquote(paste0("In ", names[1], " WGS TAT is ", round(tat_calender_this_quartile, 1), " days"))

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
invisible(dev.off())

noquote(paste0(pct[1], "% of LR patients were considered eligible for a WGS-informed treatment."))

# LR patients: who got actually treated based on WGS biomarker? ----------------------------------------------
yes_treated_based_on_biomarker <- nrow(benefit_tracking_WGS[benefit_tracking_WGS$Patient.got.treated.based.on.WGS.biomarker. == "Yes" & benefit_tracking_WGS$WGS.allowed.therapy. == "Yes" & benefit_tracking_WGS$Category == "LR", ])
no_treated_based_on_biomarker <- nrow(benefit_tracking_WGS[benefit_tracking_WGS$Patient.got.treated.based.on.WGS.biomarker. == "No" & benefit_tracking_WGS$WGS.allowed.therapy. == "Yes" & benefit_tracking_WGS$Category == "LR", ])
unknown_treated_based_on_biomarker <- nrow(benefit_tracking_WGS[benefit_tracking_WGS$Patient.got.treated.based.on.WGS.biomarker. == "Unknown" & benefit_tracking_WGS$WGS.allowed.therapy. == "Yes" & benefit_tracking_WGS$Category == "LR", ])
question_mark_treated_based_on_biomarker <- nrow(benefit_tracking_WGS[benefit_tracking_WGS$Patient.got.treated.based.on.WGS.biomarker. == "?" & benefit_tracking_WGS$WGS.allowed.therapy. == "Yes" & benefit_tracking_WGS$Category == "LR", ])
slices <- c(yes_treated_based_on_biomarker, no_treated_based_on_biomarker, unknown_treated_based_on_biomarker, question_mark_treated_based_on_biomarker)

lbls <- c("Yes:", "No:", "Unknown (definitive):", "Unknown (for now):")
pct <- round(slices/sum(slices)*100,1)
lbls <- paste(lbls, pct)
lbls <- paste0(lbls,"%")

pdf(file= paste0(wd,"treated_based_on_WGS.pdf"), width = 10, height = 7)
pie(slices,labels = lbls, col=c("lightgreen","blue","orange","red"),main="Patient got treated based on WGS biomarker?")
invisible(dev.off())

noquote(paste0(pct[1] + pct[3] + pct[4], "% of this group started (yes) or potentially started/starting (unknown) a WGS-informed treatment"))

# Impact category for LR patients NOT participating (MIGHT NEED ADJUSTMENTS IF NEW CATEGORIES ARE ADDED)  ----------------------------------------------
impact_category <- data.frame()
impact_category_columns <- c()

no_slots <- nrow(benefit_tracking_WGS[benefit_tracking_WGS$WGS.impact.category == "No slots available for WGS-informed treatment" & benefit_tracking_WGS$WGS.allowed.therapy. == "Yes" & benefit_tracking_WGS$Patient.got.treated.based.on.WGS.biomarker. == "No" & benefit_tracking_WGS$Category == "LR", ])
non_WGS_pref <- nrow(benefit_tracking_WGS[benefit_tracking_WGS$WGS.impact.category == "Non-WGS-informed treatment considered preferred" & benefit_tracking_WGS$WGS.allowed.therapy. == "Yes" & benefit_tracking_WGS$Patient.got.treated.based.on.WGS.biomarker. == "No" & benefit_tracking_WGS$Category == "LR", ])
not_fit <- nrow(benefit_tracking_WGS[benefit_tracking_WGS$WGS.impact.category == "Patient not fit enough for treatment" & benefit_tracking_WGS$WGS.allowed.therapy. == "Yes" & benefit_tracking_WGS$Patient.got.treated.based.on.WGS.biomarker. == "No" & benefit_tracking_WGS$Category == "LR", ])
not_willing <- nrow(benefit_tracking_WGS[benefit_tracking_WGS$WGS.impact.category == "Patient not willing to take treatment" & benefit_tracking_WGS$WGS.allowed.therapy. == "Yes" & benefit_tracking_WGS$Patient.got.treated.based.on.WGS.biomarker. == "No" & benefit_tracking_WGS$Category == "LR", ])
not_effective <- nrow(benefit_tracking_WGS[benefit_tracking_WGS$WGS.impact.category == "WGS-informed treatment not expected to be effective" & benefit_tracking_WGS$WGS.allowed.therapy. == "Yes" & benefit_tracking_WGS$Patient.got.treated.based.on.WGS.biomarker. == "No" & benefit_tracking_WGS$Category == "LR", ])
not_meeting_criteria <- nrow(benefit_tracking_WGS[benefit_tracking_WGS$WGS.impact.category == "Does not meet all trial criteria for WGS informed treatment" & benefit_tracking_WGS$WGS.allowed.therapy. == "Yes"  & benefit_tracking_WGS$Patient.got.treated.based.on.WGS.biomarker. == "No" & benefit_tracking_WGS$Category == "LR", ])
soc_treatment <- nrow(benefit_tracking_WGS[benefit_tracking_WGS$WGS.impact.category == "Treated using SOC for new diagnosis based on WGS" & benefit_tracking_WGS$WGS.allowed.therapy. == "Yes"  & benefit_tracking_WGS$Patient.got.treated.based.on.WGS.biomarker. == "No" & benefit_tracking_WGS$Category == "LR", ])
soc_treatment_cup <- nrow(benefit_tracking_WGS[benefit_tracking_WGS$WGS.impact.category == "Treated using SOC for unraveled diagnosis based on WGS" & benefit_tracking_WGS$WGS.allowed.therapy. == "Yes"  & benefit_tracking_WGS$Patient.got.treated.based.on.WGS.biomarker. == "No" & benefit_tracking_WGS$Category == "LR", ])

if (no_slots>0) {
  impact_category[1,'no_slots'] <- no_slots
  impact_category_columns <- append(impact_category_columns, "No slots available for WGS-informed treatment")
}

if (non_WGS_pref>0) {
  impact_category['non_WGS_pref'] <- non_WGS_pref
  impact_category_columns <- append(impact_category_columns, "Non-WGS-informed treatment considered preferred")
}

if (not_fit>0) {
  impact_category['not_fit'] <- not_fit
  impact_category_columns <- append(impact_category_columns, "Patient not fit enough for treatment")
}

if (not_willing>0) {
  impact_category['not_willing'] <- not_willing
  impact_category_columns <- append(impact_category_columns, "Patient not willing to take treatment")
}

if (not_effective>0) {
  impact_category['not_effective'] <- not_effective
  impact_category_columns <- append(impact_category_columns, "Treatment not expected to be effective")
}

if (not_meeting_criteria>0) {
  impact_category['not_meeting_criteria'] <- not_meeting_criteria
  impact_category_columns <- append(impact_category_columns, "Patient does not meet other trial criteria")
}

if (soc_treatment>0) {
  impact_category['soc_treatment'] <- soc_treatment
  impact_category_columns <- append(impact_category_columns, "Treated using SOC for new diagnosis based on WGS")
}

# Should be 0 because category belongs to CUP
if (soc_treatment_cup>0) {
  impact_category['soc_treatment_cup'] <- soc_treatment_cup
  impact_category_columns <- append(impact_category_columns, "Treated using SOC for WGS-unraveled tumor type")
  noquote("- WARN! - An invalid category has been added to the WGS impact category plot!")
}

colnames(impact_category) <- impact_category_columns

pdf(file= paste0(wd,"WGS_impact_category.pdf"), width = 20, height = 7)
par(mar=c(5,25,5,5))
barplot(as.matrix(rev(impact_category)), xlim = c(0,35), col = "blue", xlab = "Count of WGS impact category", las=1, horiz = TRUE, main = "WGS impact category if patient not treated based on WGS biomarker")
grid(nx=NULL,ny=NA,lty=1,col="gray",lwd=1)
barplot(as.matrix(rev(impact_category)), xlim = c(0,35), col = "blue", xlab = "Count of WGS impact category", las=1, horiz = TRUE, main = "WGS impact category if patient not treated based on WGS biomarker", add=TRUE)
invisible(dev.off())

lr_with_biomarker <- nrow(benefit_tracking_WGS[benefit_tracking_WGS$WGS.allowed.therapy. == "Yes" & benefit_tracking_WGS$Category == "LR", ])
lr_with_biomarker_not_treat <- nrow(benefit_tracking_WGS[benefit_tracking_WGS$WGS.allowed.therapy. == "Yes" & benefit_tracking_WGS$Category == "LR" & benefit_tracking_WGS$Patient.got.treated.based.on.WGS.biomarker. == "No", ])
lr_with_biomarker_not_treated_perc <- (lr_with_biomarker_not_treat / lr_with_biomarker) * 100
lr_with_biomarker_not_treated_not_fit_perc <- (nrow(benefit_tracking_WGS[benefit_tracking_WGS$WGS.allowed.therapy. == "Yes" & benefit_tracking_WGS$Category == "LR" & benefit_tracking_WGS$Patient.got.treated.based.on.WGS.biomarker. == "No" & benefit_tracking_WGS$WGS.impact.category == "Patient not fit enough for treatment", ]) / lr_with_biomarker_not_treat) * 100

noquote(paste0(round(lr_with_biomarker_not_treated_perc), "% of LR patients with a WGS-biomarker did not get WGS-informed treatment"))
noquote(paste0(round(lr_with_biomarker_not_treated_not_fit_perc), "% of this group of patients were not fit enough anymore for treatment"))

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
pie(slices,labels = lbls, col=c("orange","blue","red"))
invisible(dev.off())

noquote(paste0("For ", nr_confirmation_or_diagnosis, "% CUP patients WGS led to confirmation (", nr_confirmation, "%) or actual diagnosis (", nr_actual_diagnosis, "%)"))

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
invisible(dev.off())

WTS_CUP_improved_prediction_perc <- WTS_CUP_improved_prediction / nr_of_cups * 100

noquote(paste0("For ", round(WTS_CUP_improved_prediction_perc), "% of CUPs, addition of WTS yields a high confidence prediction while WGS alone does not"))

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
invisible(dev.off())

WTS_LR_potentially_actionable_perc <- WTS_LR_potential_actionable / nr_of_lr * 100

WTS_succeeded <- nrow(benefit_tracking_WTS[benefit_tracking_WTS$RNA.prep.quality.sufficient. == "Yes", ])
WTS_not_succeeded <- nrow(benefit_tracking_WTS[benefit_tracking_WTS$RNA.prep.quality.sufficient. == "No", ])
WTS_succeeded_perc <- (WTS_succeeded / (WTS_not_succeeded + WTS_succeeded))*100

WTS_hypoth_benefit_perc <- nrow(benefit_tracking_WTS[benefit_tracking_WTS$RNA.prep.quality.sufficient. == "Yes" & (benefit_tracking_WTS$WTS.detected.novel.actionable.events..not.detected.or.detectable.by.WGS. == "Yes" | benefit_tracking_WTS$WTS.improved.CUPPA.prediction.for.determining.primary.tumor.location. == "Yes" & benefit_tracking_WTS$Category == "CUP" ) , ]) / WTS_succeeded * 100

noquote(paste0("For ", round(WTS_LR_potentially_actionable_perc), "% of LRs, WTS discovers a potentially actionable biomarker not found by WGS"))
noquote(paste0("WTS succeeds for ", round(WTS_succeeded_perc), "% of biopsies, offering hypothetical benefit in ", round(WTS_hypoth_benefit_perc), "% in case of successful sequencing"))