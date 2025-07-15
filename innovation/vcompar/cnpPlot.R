#!/usr/bin/Rscript
# Author: Teoman Deger
# -----------------------./
args = commandArgs(trailingOnly=TRUE)
setwd("/home/tdeger/coloOldvNew/")

suppressMessages(library(ggplot2))
suppressMessages(library(dplyr))

#sanity checks
if (length(args)!=5) {
  stop("plotter requires 5 arguments, not all were found", call.=FALSE)
}

file1<-args[1]
file2<-args[2]
file3<-args[3]
file4<-args[4]
refGen<-args[5]

chromosome_data <- read.table("genome_length.txt", header= TRUE)
if (refGen == "hg19"){
  chromosome_data <- chromosome_data[,c("chromosome","length.hg19","length_cumsum.hg19")] 
} else if (refGen == "hg38"){
  chromosome_data <- chromosome_data[,c("chromosome","length.hg38","length_cumsum.hg38")]
} else {
  print("provide hg19 or hg38")
  q()
}
colnames(chromosome_data) <- c("chromosome", "length", "length_cumsum")

######### Read in the data
data_left <- read.table(file1, header=TRUE)
data_left <- data_left [, 1:7]
segm_left <- read.table(file2, header = TRUE)
data_right <- read.table(file3, header=TRUE)
data_right <- data_right [, 1:7]
segm_right <- read.table(file4, header = TRUE)

######### LEFT SIDE
colnames(data_left) <- c("chromosome",	"position",	"referenceReadCount",	"tumorReadCount",	"referenceGCRatio",	"tumorGCRatio", "referenceGCDiploidRatio")
data_left[(data_left$tumorGCRatio == -1),colnames(data_left)=="tumorGCRatio"] <-NA
merged_data_left <- merge(data_left, chromosome_data, by = "chromosome")
merged_data_left$adjusted_position <- merged_data_left$position + merged_data_left$length_cumsum
final_data_left <- merged_data_left[, c("chromosome", "adjusted_position", "referenceReadCount", "tumorReadCount", 
                                        "referenceGCRatio", "tumorGCRatio", "referenceGCDiploidRatio")]
colnames(segm_left) <- c("sampleID",	"chromosome",	"arm",	"start.pos",	"end.pos",	"n.probes", "mean")
merged_segm_left <- merge(segm_left, chromosome_data, by = "chromosome")
merged_segm_left$adjusted_start <- merged_segm_left$start.pos + merged_segm_left$length_cumsum
merged_segm_left$adjusted_end <- merged_segm_left$end.pos + merged_segm_left$length_cumsum
final_segm_left<- merged_segm_left[,c("chromosome", "mean", "adjusted_start", "adjusted_end")]
final_segm_left$group <- "left"

rm(data_left)
rm(merged_data_left)
rm(segm_left)
rm(merged_segm_left)

######## RIGHT SIDE
colnames(data_right) <- c("chromosome",	"position",	"referenceReadCount",	"tumorReadCount",	"referenceGCRatio",	"tumorGCRatio", "referenceGCDiploidRatio")
data_right[(data_right$tumorGCRatio == -1),colnames(data_right)=="tumorGCRatio"] <-NA
merged_data_right <- merge(data_right, chromosome_data, by = "chromosome")
merged_data_right$adjusted_position <- merged_data_right$position + merged_data_right$length_cumsum
final_data_right <- merged_data_right[, c("chromosome", "adjusted_position", "referenceReadCount", "tumorReadCount", 
                                          "referenceGCRatio", "tumorGCRatio", "referenceGCDiploidRatio")]
colnames(segm_right) <- c("sampleID",	"chromosome",	"arm",	"start.pos",	"end.pos",	"n.probes", "mean")
merged_segm_right <- merge(segm_right, chromosome_data, by = "chromosome")
merged_segm_right$adjusted_start <- merged_segm_right$start.pos + merged_segm_right$length_cumsum
merged_segm_right$adjusted_end <- merged_segm_right$end.pos + merged_segm_right$length_cumsum
final_segm_right<- merged_segm_right[,c("chromosome", "mean", "adjusted_start", "adjusted_end")]
final_segm_right$group <- "right"

rm(data_right)
rm(merged_data_right)
rm(segm_right)
rm(merged_segm_right)


### log transform
final_data_left$log2tumorGCRatio <- log2(final_data_left$tumorGCRatio)
final_data_right$log2tumorGCRatio <- log2(final_data_right$tumorGCRatio)


### Add a new column to distinguish the two dataframes
final_data_left$source <- 'dfLeft'
final_data_right$source <- 'dfRight'

### Combine the two dataframes into one
combined_df <- bind_rows(final_data_left, final_data_right)
combined_seg<- bind_rows(final_segm_left, final_segm_right)

rm(final_data_left)
rm(final_segm_left)
rm(final_data_right)
rm(final_segm_right)


### Base Dot Plot
p <- ggplot(combined_df, aes(x = adjusted_position, y = log2tumorGCRatio, color = source)) +
  geom_point(shape = 46, linewidth = 1, alpha = 0.05) + 
  labs(title = "Comparison of Outcomes",
       x = "X Axis Position",
       y = "Outcome Values") +
  theme_minimal() +
  scale_color_manual(values = c("dfLeft" = "grey69", "dfRight" = "grey32", "left" = "darkblue", "right" = "darkred")) +
  #scale_color_manual(values = c('blue', 'red'), name = "Source") +  # Custom color for points
  coord_cartesian(ylim = c(-2.5, 2.5))  # Set y-axis limits

# Add Lines for Calls
p <- p + geom_segment(data = combined_seg,
                      aes(x = adjusted_start, xend = adjusted_end, 
                          y = mean, yend = mean, color = group), 
                      size = 1.2, inherit.aes = FALSE) +  # Separate color scale for segments
  scale_color_manual(values = c("dfLeft" = "grey69", "dfRight" = "grey32", "left" = "darkblue", "right" = "darkred")) +
  #scale_color_manual(values = c("left" = "darkblue", "right" = "darkred"), name = "Group") +  # Custom color for lines
  guides(color = guide_legend(override.aes = list(linewidth = 3))) +  # Adjust legend appearance if needed
  theme(legend.position = "top") +
  labs(color = "Category")  # Legend title for group

png("CobaltDistPlot.png", height = 768, width = 1024)
plot(p)
dev.off()

q()
