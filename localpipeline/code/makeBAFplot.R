library("ggplot2")
library("gtools")
library("graphics")
library("pastecs")

args <- commandArgs()
#print(args)
outdir <- args[6]
print(outdir)
baffile <- args[7]
print(baffile)
sample <- unlist(strsplit(baffile,"/"))
print(sample)
sample <- unlist(strsplit(sample[length(sample)],"_"))[1]
print(sample)

f2 <- function(x) sum(unlist(x)[2:length(x)])
f1 <- function(x) sum(unlist(x))

determine_peaks <- function(region) {
    # ADD CHECK FOR NUMBER OF ROWS >= 50?
    quant <- quantile(region$baf)
    if(quant[2]==quant[4]) {
        # LOW COMPLEXITY
        df <- data.frame(baf=quant[2], dens=99.99)
        df$CHROM <- region$CHROM[1]
        df$POS <- mean(region$POS)
        df$CALL <- NA
        return(df)
    }

    d <- density(region$baf, bw="sj")

    ts_y<-ts(d$y)
    tp<-turnpoints(ts_y)
    peaks <- data.frame(baf=d$x[tp$tppos], dens=d$y[tp$tppos])
    df <- data.frame(subset(peaks, dens>=2.0))
    df$CHROM <- region$CHROM[1]
    df$POS <- mean(region$POS)
    df$CALL <- NA

    df <- df[order(df$dens, decreasing=T),]
    return(df)
}

rowdat <- data.frame(read.table(baffile, header=F, sep='\t'))
colnames(rowdat) <- c("CHROM","POS","CALL","baf")
rowdat$baf <- as.numeric(as.character(rowdat$baf))
rowdat$baf[is.na(rowdat$baf)] <- 0.0


#df <- subset(rowdat, CHROM %in%c(12,13))
df <- data.frame()
binsize <- 1000
chromosomes <- mixedsort(as.character(unique((rowdat$CHROM))))

for (i in c(1:length(chromosomes))) {
    chrom <- chromosomes[i]
    #print(chrom)

    tmp <- subset(rowdat, CHROM==chrom)

    for (j in seq(1, nrow(tmp), by=binsize/2)) {
        maxy <- j+binsize

        if (maxy>nrow(tmp)) {maxy<-nrow(tmp)}

        region <- tmp[j:maxy,]
        if (nrow(region) >= 100) {
            peaks <- determine_peaks(region)
            #print(nrow(peaks))
            df <- rbind(df, peaks[1:2,])
        }

        # check if HETRO / HOMOZYGOUS calls present
        if (peaks$baf[1]<=0.1 || peaks$baf[1]>=0.9) {
            # RE-RUN for intermediate BAFs
            region <- subset(region, baf>0.1 & baf<0.9)
            if (nrow(region) >= 100) {
                peaks <- determine_peaks(region)
                #print(nrow(peaks))
                df <- rbind(df, peaks)
            }
        }
    }

}

df$baf <- round(df$baf,3)
df <- df[!is.na(df$CHROM),]
df$CHROM <- factor(df$CHROM, levels=chromosomes)
#rm(rowdat)

pdf(file=paste0(outdir,sample,"_BAF.pdf"), width=15, height=2, pointsize=6, useDingbats=FALSE)

p <- ggplot(df, aes(POS, baf, group=CHROM)) + geom_point(shape=20, size=0.2, alpha=0.5) + facet_wrap(~CHROM, nrow=1, scales="free_x")

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
dev.off()
rm(p)
