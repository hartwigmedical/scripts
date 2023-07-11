#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(bizdays)) # Calculations of working days TATs

args = commandArgs(trailingOnly=TRUE)
if (length(args)==0) {
  stop("At least one argument must be supplied (tumor sample barcode).n", call.=FALSE)
}

make_business_calendar <- function() {
  non_business_weekdays <- c('saturday','sunday')
  non_business_holidays <- c("2023-01-01", "2023-04-07", "2023-04-09",
                             "2023-04-10", "2023-04-27", "2023-05-18",
                             "2023-05-29", "2023-12-25", "2023-12-26",
                             "2024-01-01", "2024-03-29", "2024-04-01",
                             "2024-04-26", "2024-05-09", "2024-05-20",
                             "2024-12-25", "2024-12-26")
  calendar <- create.calendar('nl_business_calendar',
                              weekdays=non_business_weekdays,
                              holidays=non_business_holidays)
  return(calendar)
}

barcode <- args[1]
business_calendar <- make_business_calendar()

cmd <- paste("gather_dates_for_tat.sh", barcode)
barcode_dates <- system(command=cmd, intern=T)

if (barcode_dates[2] == "NA") {
  print("Tumor material not yet received")
  q()
}
if (barcode_dates[3] == "NA") {
  print("Reference material not yet received")
  q()
}

if (as.Date(barcode_dates[2]) > as.Date(barcode_dates[3])) {
  received_date <- as.Date(barcode_dates[2]) # tumor
} else {
  received_date <- as.Date(barcode_dates[3]) # reference
}

if (barcode_dates[4] == "NA") {
  reported_date <- format(Sys.Date(), "%Y-%m-%d")
  print("This sample is still in process!")
  print("The current TAT is:")
} else {
  reported_date <- as.Date(barcode_dates[4])
  print("This OncoAct has been reported on:")
  print(reported_date)
  print("The TAT of reporting was:")
}
tat <- bizdays(received_date, reported_date, cal=business_calendar) + 1
print(tat)