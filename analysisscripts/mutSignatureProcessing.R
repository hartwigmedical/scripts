library(devtools)
library(MutationalPatterns)
library(RMySQL)
library(data.table)
library(grid)
library(gridExtra)
library("NMF")
library(ggplot2)

myCOLORS = c("#ff994b", "#463ec0", "#88c928", "#996ffb", "#68b1c0", "#e34bd9", "#106b00", "#d10073", "#98d76a",
             "#6b3a9d", "#d5c94e", "#0072e2", "#ff862c", "#31528d", "#d7003a", "#323233", "#ff4791", "#01837a",
             "#ff748a", "#777700", "#ff86be", "#4a5822", "#ffabe4", "#6a4e03", "#c6c0fb", "#ffb571", "#873659",
             "#dea185", "#a0729d", "#8a392f")

grid_arrange_shared_legend <- function(..., ncol = length(list(...)), nrow = 1, position = c("bottom", "right")) {

    plots <- list(...)
    position <- match.arg(position)
    g <- ggplotGrob(plots[[1]] + theme(legend.position = position))$grobs
    legend <- g[[which(sapply(g, function(x) x$name) == "guide-box")]]
    lheight <- sum(legend$height)
    lwidth <- sum(legend$width)
    gl <- lapply(plots, function(x) x + theme(legend.position = "none"))
    gl <- c(gl, ncol = ncol, nrow = nrow)

    combined <- switch(position,
                     "bottom" = arrangeGrob(do.call(arrangeGrob, gl),
                                            legend,
                                            ncol = 1,
                                            heights = unit.c(unit(1, "npc") - lheight, lheight)),
                     "right" = arrangeGrob(do.call(arrangeGrob, gl),
                                           legend,
                                           ncol = 2,
                                           widths = unit.c(unit(1, "npc") - lwidth, lwidth)))

    grid.newpage()
    grid.draw(combined)

    # return gtable invisibly
    invisible(combined)

}

dataFile = "~/hmf/mutSignature2.RData"
load(dataFile)

# distinct cancer types, excluding NA
cancer_types = na.omit(unique(cohort$cancerType))
cancer_types = sort(cancer_types)

# filter by cancerType
samples_with_breast_cancer = cohort[cancerType == "Breast"]$sampleId

for (s in samples_with_breast_cancer) {
    signatures[[s]]
}
