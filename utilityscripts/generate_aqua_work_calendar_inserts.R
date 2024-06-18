#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(bizdays)) # Calculations of working days TATs

args <- commandArgs(trailingOnly=TRUE)
if (length(args)!=4) {
  stop("Supply four arguments: start_date, end_date, holiday_file, output_file", call.=FALSE)
}

start_date <- args[1]
end_date <- args[2]
holiday_file <- args[3]
output_file <- args[4]

# Check start_date
check_start_date <- try(as.Date(start_date, format="%Y-%m-%d"))
if("try-error" %in% class(check_start_date) || is.na(check_start_date)) {
  print("ERROR: start_date is not a valid date")
  quit()
}

# Check end_date
check_end_date <- try(as.Date(end_date, format="%Y-%m-%d"))
if("try-error" %in% class(check_end_date) || is.na(check_end_date)) {
  print("ERROR: end_date is not a valid date")
  quit()
}

# Check start_date greater than end_date for validity
if (as.Date(start_date, format="%Y-%m-%d") >= as.Date(end_date, format="%Y-%m-%d")) {
  print("ERROR: start_date is greater than end_date")
  quit()
}

# Check input holiday_file
if (! file.exists(holiday_file)) {
  print("ERROR: cannot find holiday_file")
  quit()
}

# Load holidays and make calendar
weekend_list <- c("saturday", "sunday")
holidays_list <- scan(holiday_file, character())
work_calendar <- create.calendar('nl_business_calendar',
                                 weekdays=weekend_list,
                                 holidays=holidays_list,
                                 start.date=start_date,
                                 end.date=end_date
)

# Loop through all dates and determine if it is a business day, write the combination to output file
fileConn <- file(output_file, "w")
writeLines(c("INSERT INTO work_calender (date_time, is_working_day) VALUES"), fileConn)
for (date in seq(as.Date(start_date), as.Date(end_date)-1, by="day")) {
  writeLines(paste0( "('", as.Date(date), "',", as.integer(is.bizday(date, work_calendar)), ")," ), fileConn)
}
writeLines(paste0( "('", as.Date(end_date), "',", as.integer(is.bizday(end_date, work_calendar)), ");" ), fileConn)
close(fileConn)