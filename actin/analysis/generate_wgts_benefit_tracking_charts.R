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
Q1_2024 <- c(240101, 240331)
Q2_2024 <- c(240401, 240630)
Q3_2024 <- c(240701, 240930)
Q4_2024 <- c(241001, 241231)

#add new quartile
names <- c("2024 Q4", "2024 Q3", "2024 Q2", "2024 Q1", "2023 Q4", "2023 Q3", "2023 Q2", "2023 Q1", "2022 Q4", "2022 Q3", "2022 Q2", "2022 Q1", "2021 Q4", "2021 Q3")
quartiles <- list(Q4_2024, Q3_2024, Q2_2024, Q1_2024, Q4_2023, Q3_2023, Q2_2023, Q1_2023, Q4_2022, Q3_2022, Q2_2022, Q1_2022, Q4_2021, Q3_2021)


# Getting started ---------------------------------------------------------
setwd(wd)
filename_WGS <- "EMC ACTIN WGS+WTS Benefit Tracking - WGS.tsv"
filename_WTS <- "EMC ACTIN WGS+WTS Benefit Tracking - WTS.tsv"
benefit_tracking_WGS <- read.table(file=filename_WGS, sep = '\t', header = TRUE)
benefit_tracking_WTS <- read.table(file=filename_WTS, sep = '\t', header = TRUE)


# LR/CUP patients: Nr of WGS analyses  ------------------------------------------------------------------
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
total_nr_patients <- subset(benefit_tracking_WGS, benefit_tracking_WGS$WGS.report.date <= quartiles[[1]][2] & benefit_tracking_WGS$WGS.report.date >= quartiles[[length(quartiles)]][1] & !is.nan(benefit_tracking_WGS$WGS.report.date), select = Patient.ID)
total_nr_patients_value <- length(unique(total_nr_patients$Patient.ID))

pdf(file= paste0(wd,"LR_CUP.pdf"), width = 10, height = 7)
par(mar=c(5,5,5,5),mfrow=c(1,1))
barplot(as.matrix(LR_CUP), xlim = c(0,max(LR_CUP[1,])+max(LR_CUP[2,])+5), col=c("blue","red"), las=1, legend = c("LR", "CUP"), horiz = TRUE, args.legend = list(x ='topright', inset = c(-0.05,0.05)))
grid(nx=NULL,ny=NA,lty=1,col="gray",lwd=1)
barplot(as.matrix(LR_CUP), xlim = c(0,max(LR_CUP[1,])+max(LR_CUP[2,])+5), col=c("blue","red"), las=1, legend = c("LR", "CUP"), horiz = TRUE, args.legend = list(x ='topright', inset = c(-0.05,0.05)), add=TRUE)
invisible(dev.off())

noquote("LR_CUP.pdf")
noquote(paste0("Overall, ", total_nr_biopsies, " biopsies (", total_nr_biopsies_LR, " LR, ", total_nr_biopsies_CUP, " CUP) for ", total_nr_patients_value, " patients have been evaluated in ACTIN "))
noquote(paste0("In ", names[1], ", ", LR_CUP[1,1] + LR_CUP[2,1], " WGS analyses were reported (", LR_CUP[1,1], " LR and ", LR_CUP[2,1], " CUP)"))
noquote("")

# LR & CUP patients: Sufficient tumor cells --------------------------------------------------
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

sufficient_nr_total <- nrow(benefit_tracking_WGS[benefit_tracking_WGS$WGS.performed.successfully. == "Yes", ])
sufficient_pct_total <- sufficient_nr_total/total_nr_biopsies*100

pdf(file= paste0(wd,"sufficient_tumor_cells.pdf"), width = 10, height = 7)
par(mar=c(5,5,5,5)+.1)
barplot(as.matrix(suff_cells), xlim = c(0,100), col=c("blue","red"), las=1, horiz = TRUE, xaxt = "n")
grid(nx=NULL,ny=NA,lty=1,col="gray",lwd=1)
barplot(as.matrix(suff_cells), xlim = c(0,100), col=c("blue","red"), las=1, horiz = TRUE, add=TRUE, xaxt = "n")
axis(1,at=c(0,20,40,60,80,100),labels=paste0(c(0,20,40,60,80,100), "%"))
invisible(dev.off())

noquote("sufficient_tumor_cells.pdf")
noquote(paste0("Overall, ", round(sufficient_pct_total), "% of biopsies contain sufficient tumor cells for WGTS analysis"))
noquote(paste0("In ", names[1], ", ", round(suff_cells[1,1]), "% of biopsies contained sufficient tumor cells"))
noquote("")

# LR patients: Sufficient tumor cells --------------------------------------------------
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

sufficient_nr_LR_total <- nrow(benefit_tracking_WGS[benefit_tracking_WGS$Category == "LR" & benefit_tracking_WGS$WGS.performed.successfully. == "Yes", ])
sufficient_pct_LR_total <- sufficient_nr_LR_total/total_nr_biopsies_LR*100

pdf(file= paste0(wd,"sufficient_tumor_cells_LR.pdf"), width = 10, height = 7)
par(mar=c(5,5,5,5)+.1)
barplot(as.matrix(suff_cells_LR,1), xlim = c(0,100), col=c("blue","red"), las=1, horiz = TRUE, main = "LR", xaxt = "n")
grid(nx=NULL,ny=NA,lty=1,col="gray",lwd=1)
barplot(as.matrix(suff_cells_LR,1), xlim = c(0,100), col=c("blue","red"), las=1, horiz = TRUE, main = "LR", add=TRUE, xaxt = "n")
axis(1,at=c(0,20,40,60,80,100),labels=paste0(c(0,20,40,60,80,100), "%"))
invisible(dev.off())

noquote("sufficient_tumor_cells_LR.pdf")
noquote(paste0("Overall, ", round(sufficient_pct_LR_total), "% of biopsies contain sufficient tumor cells for WGS analysis for LR"))
noquote(paste0("In ", names[1], ", " , round(suff_cells_LR[1,1]), "% of biopsies contain sufficient tumor cells for WGS analysis for LR"))
noquote("")

# CUP patients: Sufficient tumor cells --------------------------------------------------
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

sufficient_nr_CUP_total <- nrow(benefit_tracking_WGS[benefit_tracking_WGS$Category == "CUP" & benefit_tracking_WGS$WGS.performed.successfully. == "Yes", ])
sufficient_pct_CUP_total <- sufficient_nr_CUP_total/total_nr_biopsies_CUP*100

pdf(file= paste0(wd,"sufficient_tumor_cells_CUP.pdf"), width = 10, height = 7)
par(mar=c(5,5,5,5)+.1)
barplot(as.matrix(suff_cells_CUP,1), xlim = c(0,100), col=c("blue","red"), las=1, horiz = TRUE, main = "CUP", xaxt = "n")
grid(nx=NULL,ny=NA,lty=1,col="gray",lwd=1)
barplot(as.matrix(suff_cells_CUP,1), xlim = c(0,100), col=c("blue","red"), las=1, horiz = TRUE, main = "CUP", add=TRUE, xaxt = "n")
axis(1,at=c(0,20,40,60,80,100),labels=paste0(c(0,20,40,60,80,100), "%"))
invisible(dev.off())

noquote("sufficient_tumor_cells_CUP.pdf")
noquote(paste0("Overall, ", round(sufficient_pct_CUP_total), "% of biopsies contain sufficient tumor cells for WGS analysis for CUP"))
noquote(paste0("In ", names[1], ", " ,round(suff_cells_CUP[1,1]), "% of biopsies contain sufficient tumor cells for WGS analysis for CUP"))
noquote("")

# TAT in calendar days ---------------------------------------------------------------------
tat_calendar <- data.frame()
i <- 1

for (x in quartiles) {
  tat <- subset(benefit_tracking_WGS, benefit_tracking_WGS$WGS.report.date <= x[2] & benefit_tracking_WGS$WGS.report.date >= x[1] & !is.nan(benefit_tracking_WGS$WGS.report.date), select=TAT...calendar.days.)
  tat_calendar[1,i] <- mean(tat$TAT...calendar.days.)
  i <- i+1
}
colnames(tat_calendar) <- names

average_tat <- mean(benefit_tracking_WGS$TAT...calendar.days.)

pdf(file= paste0(wd,"TAT_calendar_days.pdf"), width = 10, height = 7)
barplot(as.matrix(rev(tat_calendar)), ylab = "TAT calendar days", ylim = c(0,max(tat_calendar)+1), col = "blue", las = 2)
grid(nx=NA,ny=NULL,lty=1,col="gray",lwd=1)
barplot(as.matrix(rev(tat_calendar)), ylab = "TAT calendar days", ylim = c(0,max(tat_calendar)+1), col = "blue", add=TRUE, las = 2)
invisible(dev.off())

tat_calendar_this_quartile <- tat_calendar[1, 1]

noquote("TAT_calendar_days.pdf")
noquote(paste0("Average TAT (turnaround-time) in calendar days for WGS is ", round(average_tat,1), " days"))
noquote(paste0("In ", names[1], ", WGS average TAT is ", round(tat_calendar_this_quartile, 1), " days"))
noquote("")

# LR patients: WGS allowed therapy? ----------------------------------------------
allowed_therapy <- nrow(benefit_tracking_WGS[benefit_tracking_WGS$WGS.allowed.therapy. == "Yes" & benefit_tracking_WGS$Category == "LR" , ])
non_allowed_therapy <- nrow(benefit_tracking_WGS[benefit_tracking_WGS$WGS.allowed.therapy. == "No" & benefit_tracking_WGS$Category == "LR", ])
slices <- c(allowed_therapy, non_allowed_therapy)

lbls <- c("Yes:", "No:")
pct <- round(slices/sum(slices)*100,1)
lbls <- paste(lbls, pct)
lbls <- paste0(lbls,"%")

pdf(file= paste0(wd,"WGS_allowed_therapy.pdf"), width = 10, height = 7)
pie(slices, labels = lbls, col=c("green4","red"),main="WGS allowed therapy?")
invisible(dev.off())

noquote("WGS_allowed_therapy.pdf")
noquote(paste0(pct[1], "% of LR patients with successful WGS were considered eligible for a WGS-informed treatment."))
noquote("")

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
pie(slices,labels = lbls, col=c("green4","red","darkgrey","lightgrey"),main="Patient got treated based on WGS biomarker?")
invisible(dev.off())

noquote("treated_based_on_WGS.pdf")
noquote(paste0(round(pct[1] + pct[3] + pct[4]), "% of this group actually started (yes) or potentially started/starting (unknown) a WGS-informed treatment"))
noquote("")

# LR patients: Impact category for patients NOT participating (MIGHT NEED ADJUSTMENTS IF NEW CATEGORIES ARE ADDED)  ----------------------------------------------
impact_category <- data.frame()
impact_category_columns <- c()

no_slots <- nrow(benefit_tracking_WGS[benefit_tracking_WGS$WGS.impact.category == "No slots available for WGS-informed treatment" & benefit_tracking_WGS$WGS.allowed.therapy. == "Yes" & benefit_tracking_WGS$Patient.got.treated.based.on.WGS.biomarker. == "No" & benefit_tracking_WGS$Category == "LR", ])
non_WGS_pref <- nrow(benefit_tracking_WGS[benefit_tracking_WGS$WGS.impact.category == "Non-WGS-informed treatment considered preferred" & benefit_tracking_WGS$WGS.allowed.therapy. == "Yes" & benefit_tracking_WGS$Patient.got.treated.based.on.WGS.biomarker. == "No" & benefit_tracking_WGS$Category == "LR", ])
not_fit <- nrow(benefit_tracking_WGS[benefit_tracking_WGS$WGS.impact.category == "Patient not fit enough for treatment" & benefit_tracking_WGS$WGS.allowed.therapy. == "Yes" & benefit_tracking_WGS$Patient.got.treated.based.on.WGS.biomarker. == "No" & benefit_tracking_WGS$Category == "LR", ])
not_willing <- nrow(benefit_tracking_WGS[benefit_tracking_WGS$WGS.impact.category == "Patient not willing to take treatment" & benefit_tracking_WGS$WGS.allowed.therapy. == "Yes" & benefit_tracking_WGS$Patient.got.treated.based.on.WGS.biomarker. == "No" & benefit_tracking_WGS$Category == "LR", ])
not_effective <- nrow(benefit_tracking_WGS[benefit_tracking_WGS$WGS.impact.category == "WGS-informed treatment not expected to be effective" & benefit_tracking_WGS$WGS.allowed.therapy. == "Yes" & benefit_tracking_WGS$Patient.got.treated.based.on.WGS.biomarker. == "No" & benefit_tracking_WGS$Category == "LR", ])
not_meeting_criteria <- nrow(benefit_tracking_WGS[benefit_tracking_WGS$WGS.impact.category == "Does not meet all trial criteria for WGS informed treatment" & benefit_tracking_WGS$WGS.allowed.therapy. == "Yes"  & benefit_tracking_WGS$Patient.got.treated.based.on.WGS.biomarker. == "No" & benefit_tracking_WGS$Category == "LR", ])
soc_treatment <- nrow(benefit_tracking_WGS[benefit_tracking_WGS$WGS.impact.category == "Treated using SOC for new diagnosis based on WGS" & benefit_tracking_WGS$WGS.allowed.therapy. == "Yes"  & benefit_tracking_WGS$Patient.got.treated.based.on.WGS.biomarker. == "No" & benefit_tracking_WGS$Category == "LR", ])

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

# Check if category adjustments are needed
if (sum(no_slots, non_WGS_pref, not_fit, not_willing, not_effective, not_meeting_criteria, soc_treatment) != no_treated_based_on_biomarker) {
  noquote ("- !WARN! - Not all categories could be assigned to the WGS impact category plots!")
} else {
  noquote ("- INFO - All categories are assigned, no need to add new categories to WGS_impact_category.pdf")
}

colnames(impact_category) <- impact_category_columns
order_indices <- order(as.numeric(impact_category), decreasing=FALSE)
impact_category_ordered <- impact_category[, order_indices]

pdf(file= paste0(wd,"WGS_impact_category.pdf"), width = 20, height = 7)
par(mar=c(5,25,5,5))
barplot(as.matrix(impact_category_ordered), xlim = c(0,35), col = "blue", xlab = "Number of patients", las=1, horiz = TRUE, main = "Reason if patient not treated based on WGS biomarker")
grid(nx=NULL,ny=NA,lty=1,col="gray",lwd=1)
barplot(as.matrix(impact_category_ordered), xlim = c(0,35), col = "blue", xlab = "Number of patients", las=1, horiz = TRUE, main = "Reason if patient not treated based on WGS biomarker", add=TRUE)
invisible(dev.off())

lr_with_biomarker <- nrow(benefit_tracking_WGS[benefit_tracking_WGS$WGS.allowed.therapy. == "Yes" & benefit_tracking_WGS$Category == "LR", ])
lr_with_biomarker_not_treated <- nrow(benefit_tracking_WGS[benefit_tracking_WGS$WGS.allowed.therapy. == "Yes" & benefit_tracking_WGS$Category == "LR" & benefit_tracking_WGS$Patient.got.treated.based.on.WGS.biomarker. == "No", ])
lr_with_biomarker_not_treated_perc <- (lr_with_biomarker_not_treated / lr_with_biomarker) * 100
lr_with_biomarker_not_treated_not_fit_perc <- (nrow(benefit_tracking_WGS[benefit_tracking_WGS$WGS.allowed.therapy. == "Yes" & benefit_tracking_WGS$Category == "LR" & benefit_tracking_WGS$Patient.got.treated.based.on.WGS.biomarker. == "No" & benefit_tracking_WGS$WGS.impact.category == "Patient not fit enough for treatment", ]) / lr_with_biomarker_not_treated) * 100

noquote("WGS_impact_category.pdf")
noquote(paste0(round(lr_with_biomarker_not_treated_perc), "% of LR patients with a WGS-biomarker did not get WGS-informed treatment"))
noquote(paste0(round(lr_with_biomarker_not_treated_not_fit_perc), "% of this group of patients were not fit enough anymore for treatment"))
noquote("")

# CUP patients: WGS ----------------------------------------------
cup_biopsy_contained_insuff_cells <- nrow(benefit_tracking_WGS[benefit_tracking_WGS$WGS.performed.successfully. == "No" & benefit_tracking_WGS$Category == "CUP", ])
cup_no_impact_diagnosis <- nrow(benefit_tracking_WGS[benefit_tracking_WGS$WGS.led.to..change.of..diagnosis. == "No" & benefit_tracking_WGS$WGS.led.to.confirmation.of.uncertain.diagnosis. != "Yes" & benefit_tracking_WGS$Category == "CUP", ])
cup_confirmation_diagnosis <- nrow(benefit_tracking_WGS[benefit_tracking_WGS$WGS.led.to..change.of..diagnosis. == "No" & benefit_tracking_WGS$WGS.led.to.confirmation.of.uncertain.diagnosis. == "Yes" & benefit_tracking_WGS$Category == "CUP", ])
cup_actual_diagnosis <- nrow(benefit_tracking_WGS[benefit_tracking_WGS$WGS.led.to..change.of..diagnosis. == "Yes" & benefit_tracking_WGS$Category == "CUP", ])
cup_unknown_impact_diagnosis <- nrow(benefit_tracking_WGS[benefit_tracking_WGS$WGS.led.to..change.of..diagnosis. == "Unknown" & benefit_tracking_WGS$Category == "CUP", ])

if (sum(cup_biopsy_contained_insuff_cells, cup_no_impact_diagnosis, cup_confirmation_diagnosis, cup_unknown_impact_diagnosis, cup_actual_diagnosis) != total_nr_biopsies_CUP) {
  noquote ("- !WARN! - Not all patients could be assigned to the CUP diagnosis categories!")
} else {
  noquote ("- INFO - All patients could be assigned to the CUP diagnosis categories")
}

slices <- c(cup_no_impact_diagnosis, cup_biopsy_contained_insuff_cells, cup_confirmation_diagnosis, cup_actual_diagnosis, cup_unknown_impact_diagnosis)
lbls <- c("WGS did not impact diagnosis:", "Biopsy contained insufficient tumor cells:", "WGS led to confirmation of diagnosis:", "WGS led to actual diagnosis:", "Unknown:")
pct <- round(slices/sum(slices)*100,1)
lbls <- paste(lbls, round(pct))
lbls <- paste0(lbls, "%")

pdf(file= paste0(wd,"WGS_CUP.pdf"), width = 10, height = 7)
par(mar=c(5,7,5,5))
pie(slices,labels = lbls, col=c("red","blue","lightgreen", "green4", "grey"))
invisible(dev.off())

pct_confirmation_or_diagnosis <- sum(round(pct[3]), round(pct[4]))

noquote("WGS_CUP.pdf")
noquote(paste0("For ", pct_confirmation_or_diagnosis, "% of all CUP patients WGS led to actual diagnosis or confirmation of anticipated diagnosis"))
noquote("")

# LR patients: WTS ----------------------------------------------
wts_lr_sufficient_rna <- nrow(benefit_tracking_WTS[benefit_tracking_WTS$RNA.prep.quality.sufficient. == "Yes" & benefit_tracking_WTS$Category == "LR", ])
wts_lr_potentially_actionable_and_sufficient_rna <- nrow(benefit_tracking_WTS[benefit_tracking_WTS$RNA.prep.quality.sufficient. == "Yes" & benefit_tracking_WTS$WTS.detected.novel.actionable.events..not.detected.or.detectable.by.WGS. == "Yes" & benefit_tracking_WTS$Category == "LR", ])

WTS_LR <- data.frame()
WTS_LR[1,1] <- wts_lr_potentially_actionable_and_sufficient_rna / wts_lr_sufficient_rna * 100
WTS_LR[2,1] <- 100 - wts_lr_potentially_actionable_and_sufficient_rna / wts_lr_sufficient_rna * 100
colnames(WTS_LR) <- ""

pdf(file=paste0(wd,"WTS_LR.pdf"), width = 10, height = 3)
barplot(as.matrix(WTS_LR), col=c("blue","red"), las=1, horiz = TRUE, legend = c("Yes", "No"), args.legend = list(x ='right', inset = c(1,1)), main = "LR patients with potential additional biomarker from WTS", xaxt = "n")
grid(nx=NULL,ny=NA,lty=1,col="gray",lwd=1)
barplot(as.matrix(WTS_LR), col=c("blue","red"), las=1, horiz = TRUE, legend = c("Yes", "No"), args.legend = list(x ='right', inset = c(1,1)), main = "LR patients with potential additional biomarker from WTS", add=TRUE, xaxt = "n")
axis(1,at=c(0,20,40,60,80,100),labels=paste0(c(0,20,40,60,80,100), "%"))
invisible(dev.off())

pct_wts_lr_potentially_actionable_and_sufficient_rna <- wts_lr_potentially_actionable_and_sufficient_rna / wts_lr_sufficient_rna * 100

noquote("WTS_LR.pdf")
noquote(paste0("For ", round(pct_wts_lr_potentially_actionable_and_sufficient_rna), "% of LRs, WTS discovers a potentially actionable biomarker not found or detectable by WGS"))
noquote("")

# CUP patients: WTS ----------------------------------------------
wts_cup_sufficient_rna <- nrow(benefit_tracking_WTS[benefit_tracking_WTS$RNA.prep.quality.sufficient. == "Yes" & benefit_tracking_WTS$Category == "CUP", ])
wts_cup_improved_prediction_and_sufficient_rna <- nrow(benefit_tracking_WTS[benefit_tracking_WTS$RNA.prep.quality.sufficient. == "Yes" & benefit_tracking_WTS$WTS.improved.CUPPA.prediction.for.determining.primary.tumor.location. == "Yes" & benefit_tracking_WTS$Category == "CUP", ])

WTS_CUP <- data.frame()
WTS_CUP[1,1] <- wts_cup_improved_prediction_and_sufficient_rna / wts_cup_sufficient_rna * 100
WTS_CUP[2,1] <- 100 - wts_cup_improved_prediction_and_sufficient_rna / wts_cup_sufficient_rna * 100
colnames(WTS_CUP) <- ""

pdf(file=paste0(wd,"WTS_CUP.pdf"), width = 10, height = 3)
barplot(as.matrix(WTS_CUP), col=c("blue","red"), las=1, horiz = TRUE, legend = c("Yes", "No"), args.legend = list(x ='right', inset = c(1,1)), xaxt = "n", main = "CUPs where WTS leads to high confidence CUPPA prediction")
grid(nx=NULL,ny=NA,lty=1,col="gray",lwd=1)
barplot(as.matrix(WTS_CUP), col=c("blue","red"), las=1, horiz = TRUE, legend = c("Yes", "No"), args.legend = list(x ='right', inset = c(1,1)), main = "CUPs where WTS leads to high confidence CUPPA prediction", add=TRUE, xaxt = "n")
axis(1,at=c(0,20,40,60,80,100),labels=paste0(c(0,20,40,60,80,100), "%"))
invisible(dev.off())

pct_wts_cup_improved_prediction_and_sufficient_rna <- wts_cup_improved_prediction_and_sufficient_rna / wts_cup_sufficient_rna * 100

noquote("WTS_CUP.pdf")
noquote(paste0("For ", round(pct_wts_cup_improved_prediction_and_sufficient_rna), "% of CUPs, addition of WTS yields a high confidence tumor location prediction while WGS alone does not"))
noquote("")

# LR+CUP patients: WTS added benefit  ----------------------------------------------
wts_succeeded <- nrow(benefit_tracking_WTS[benefit_tracking_WTS$RNA.prep.quality.sufficient. == "Yes", ])
wts_not_succeeded <- nrow(benefit_tracking_WTS[benefit_tracking_WTS$RNA.prep.quality.sufficient. == "No", ])
wts_succeeded_perc <- (wts_succeeded / (wts_not_succeeded + wts_succeeded))*100

wts_hypoth_benefit_perc <- nrow(benefit_tracking_WTS[benefit_tracking_WTS$RNA.prep.quality.sufficient. == "Yes" & (benefit_tracking_WTS$WTS.detected.novel.actionable.events..not.detected.or.detectable.by.WGS. == "Yes" | (benefit_tracking_WTS$WTS.improved.CUPPA.prediction.for.determining.primary.tumor.location. == "Yes" & benefit_tracking_WTS$Category == "CUP" )), ]) / wts_succeeded * 100

noquote("WTS succeeds & added benefit message:")
noquote(paste0("WTS succeeds for ", round(wts_succeeded_perc), "% of biopsies, offering hypothetical benefit (new biomarker or improved CUPPA prediction for CUPs) in ", round(wts_hypoth_benefit_perc), "% for WGS/LR patients in case of successful sequencing"))
noquote("")
noquote("Done!")