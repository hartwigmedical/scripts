library("ggplot2")
library("gtools")

myCols <- c("darkgray","steelblue3","red3","lightgray")
names(myCols) <- c("neutral1","gain","loss","neutral2")


# cat makeKaryotype.R | R --slave --args 2 24 4 500000 ".$sample_bam_name."_ratio.txt\n";
args <- commandArgs()

# Settings
ploidy <- type.convert(args[4])
nrChromsToplot <- type.convert(args[5])
maxLevelToPlot <- type.convert(args[6])
binsize <- type.convert(args[7])

# File to convert
ratiofile <- args[8]


make_plot <- function(filename) {
	df <- data.frame()

	ratio<-data.frame(read.table(filename, sep='\t', header=TRUE))

	# Clean up ratio data
	ratio <- subset(ratio, Ratio!=-1)
	ratio <- ratio[!is.na(ratio$Ratio),]

	chromosomes <- mixedsort(as.character(unique((ratio$Chromosome))))
	for (i in c(1:nrChromsToplot)) {
		chrom <- chromosomes[i]

		# Select chromosome
		tmp <- subset( ratio, Chromosome==chrom)

		# Binning of positions
		grouping <- cut(tmp$Start, c(seq(0,max(tmp$Start),binsize), max(tmp$Start)), labels=seq(binsize/2, (max(tmp$Start)+(binsize/2)), binsize))

		# Calculate median ratio
		medians <- data.frame(tapply(tmp$Ratio, grouping, median))
		rm(tmp)

		colnames(medians) <- "Ratio"
		# Alternate chromosome banding
		if ((i%%2) == 0) {
			medians$Cols <- 'neutral1'
		} else {
			medians$Cols <- 'neutral2'
		}

		# Create clean data farme
		medians <- medians[!is.na(medians$Ratio),]
		medians$Chromosome <- chrom
		medians$GenomicPosition <- as.numeric(as.character(rownames(medians)))
		medians$CopyNumber <- medians$Ratio*ploidy

		# Add event coloring
		medians$Cols[medians$CopyNumber>=ploidy+0.8] <- "gain"
		medians$Cols[medians$CopyNumber<=ploidy-0.8] <- "loss"

		# Merge data
		df <- rbind(df, medians)
	}

	df$CopyNumber[df$CopyNumber > maxLevelToPlot] <- maxLevelToPlot
	df$Chromosome <- factor(df$Chromosome, levels=chromosomes)
	rm(ratio)

	#
	# PLot the data
	pdf(file=gsub(".txt","_karyotype.pdf",filename), width=10, height=2, pointsize=6, useDingbats=FALSE)

  	p <- ggplot(df, aes(GenomicPosition, CopyNumber, group=Chromosome)) + geom_point(aes(colour=Cols), shape=20, size=0.2, alpha=0.5) + ylim(0,maxLevelToPlot) + facet_wrap(~Chromosome, nrow=1, scales="free_x") + scale_colour_manual(name="event", values=myCols)

  	# Make it look clean
  	p <- p + theme(
  			  axis.line.x=element_blank(),
  	          axis.text.x=element_blank(),
  	          axis.ticks.x=element_blank(),

  	          panel.grid.minor=element_blank(),
  	          panel.grid.major.x=element_blank(),
  	          panel.grid.major.y=element_line(colour="grey80", linetype=2, size=0.2),

  	          legend.position="none",
  	          panel.background=element_blank(),
  	          panel.border=element_blank(),
  	          plot.background=element_blank())

  	print(p)
	rm(df)

	dev.off()
}

# Make plot
make_plot(ratiofile)

