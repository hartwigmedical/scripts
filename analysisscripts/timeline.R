# This script needs a concatenated file of job logging, which should be tab-separated:
#
# (Start|End)\tJob\tTime\tStep\tHost
#
# Something like this should work:
#
# find $RUN_DIR -type f -name *.log -not -name plink.log \
#     | xargs grep -l ^Start \
#     | xargs grep -hE '^(Start|End)' \
#     > /tmp/timeline.log
#
# Some corrupt lines appear in the files. The above should filter them out, resulting in
# some jobs missing start and end times and thus extending over the whole duration. This
# would not happen on a local filesystem (you can write concurrently to the same file from
# different processes if the write is small) but Schuberg put output on NFS.
#
# At *least* until v1.12 there is not a consistent tab-separated format and manual fix-up
# is required. The script fixes the missing tab between Start/End and the job name but
# nothing else. After v1.12 remove the marked line below and add a column.
#
# Jobs need to use the same step and host to match up start and end properly. Also not the
# case for all jobs :( Could discard host. Step names should be fixed after v1.12.
#
# It is obviously quite easy to run this script from a more reliable source by taking
# database logs or any other centralised format and writing it to a file in the same
# format as expected here.
#
# Time-zone and strptime format are hard-coded for now, need altering when format changes
# for different locales. lubridate package can help with this.

library(magrittr)
library(tidyr)
library(dplyr)
library(RColorBrewer)
library(plotly)
library(htmlwidgets)
library(stringr)

getPalette <- function(categories) {
    # maximum colors of a certain type
    # num_required_colors <- length(unique(categories))
    # qualitative_palettes = brewer.pal.info[brewer.pal.info$category == 'qual', ]
    # palette <- unlist(mapply(brewer.pal, qualitative_palettes$maxcolors, rownames(qualitative_palettes)))
    # palette[1:num_required_colors]

    # pick a set and cycle
    palette_name <- "Set3"
    brewer.pal(brewer.pal.info[palette_name, "maxcolors"], name = palette_name)
}

makePlot <- function(df) {
    palette <- getPalette(df$job)
    p <- plot_ly(type = "scatter", mode = "markers")
    for (i in 1:nrow(df)) {
        p <- add_trace(p,
                       x = c(if (is.na(df$Start[i])) min(df$Start, na.rm = TRUE) else df$Start[i],
                             if (is.na(df$End[i])) max(df$End, na.rm = TRUE) else df$End[i]),
                       y = c(i, i),
                       mode = "lines",
                       line = list(color = palette[as.integer(df$job[i]) %% length(palette) + 1], width = 5),
                       showlegend = F,
                       hoverinfo = "text",
                       text = paste("Job: ", df$job[i], "<br>",
                                    "Step: ", df$step[i], "<br>",
                                    "Duration: ", round(difftime(df$End[i], df$Start[i], units = "mins"), digits = 1), " minutes<br>",
                                    "Host: ", df$host[i], "<br>")
        )
    }

    p <- p %>% layout(
        xaxis = list(showgrid = F, tickfont = list(color = "#e6e6e6")),
        annotations = list(
            list(xref = "paper", yref = "paper",
                 x = 0.85, y = 0.2,
                 text = paste0("Total Duration: ", round(difftime(max(df$End, na.rm = TRUE), min(df$Start, na.rm = TRUE), units = "hours"), digits = 1), " hours<br>",
                               "Total Time: ", round(sum(difftime(df$End, df$Start, units = "hours"), na.rm = TRUE), digits = 1), " hours<br>",
                               "Total Operations: ", length(unique(df$job)), "<br>",
                               "Total Jobs: ", nrow(df), "<br>"),
                 font = list(color = "#ffffff", size = 12),
                 ax = 0, ay = 0,
                 align = "left")),
        plot_bgcolor = "#333333",
        paper_bgcolor = "#333333"
    )
}

setClass("myDate")
setAs("character", "myDate", function(from) as.POSIXct(str_c(str_sub(from, 0, -9),
                                                             str_sub(from, -4, -1)),
                                                       format = "%a %b %e %H:%M:%S %Y",
                                                       tz = "CET"))
df <- read.csv("/tmp/timeline-clean.log",
               sep = "\t",
               header = FALSE,
               # enable after v1.12
               # col.names = c("event", "job", "time", "step", "host"),
               # colClasses = c("factor", "factor", "myDate", "factor", "factor")) %>%
               # enable after v1.12
               # remove this after v1.12
               col.names = c("event", "time", "step", "host"),
               colClasses = c("factor", "myDate", "factor", "factor")) %>%
    separate(event, c("event", "job"), " ", extra = "merge") %>%
    mutate(job = as.factor(job)) %>%
    # remove this after v1.12
    spread(event, time) %>%
    arrange(Start)

p <- makePlot(df)
htmlwidgets::saveWidget(p, "timeline.html")
p
