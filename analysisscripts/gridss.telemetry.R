library(tidyverse)
library(rtracklayer)
library(GenomeInfoDb)
library(BSgenome)
genomehg19 <- getBSgenome("hg19")

dir <- "D://hartwig//events/"
fdf <- bind_rows(lapply(list.files(path=dir, pattern="*.events.csv"), function(filename) {
	read.csv(paste0(dir, filename), header=FALSE)
}))
#fdf = read.csv("D://hartwig//CPCT02100145.gridss.assembly.bam.events.csv", header=FALSE)
names(fdf) = c("chunk", "direction", "op", "seqname", "start", "end", "records", "filtered", "usec")
ldf = fdf %>% filter(op=="load")

ggplot(ldf) +
	aes(x=usec) +
	scale_x_log10() +
	scale_y_log10() +
	geom_histogram(bins=100)

ggplot(ldf %>%
			 	mutate(bin=round(start, -5)) %>%
			 	group_by(seqname, bin) %>%
			 	summarise(
			 		usec=max(usec),
			 		records=max(records)
			 	)) +
	aes(x=bin, y=usec) +
	geom_point() +
	scale_y_log10() +
	facet_wrap(~ seqname)




binsize = 1000
# binned timing
ldfbin = ldf %>%
	mutate(start=floor(start/ binsize) * binsize) %>%
	group_by(seqname, start) %>%
	summarise(
		total_usec=sum(usec),
		max_records=max(records))
lbingr = GRanges(seqnames=ldfbin$seqname, ranges=IRanges(start=ldfbin$start, width=binsize), score=ldfbin$total_usec, total_usec=ldfbin$total_usec, max_records=ldfbin$max_records)
genome(lbingr) = "hg19"
seqlengths(lbingr) <- c(249250621,
                        135534747,
                        135006516,
                        133851895,
                        115169878,
                        107349540,
                        102531392,
                        90354753,
                        81195210,
                        78077248,
                        59128983,
                        243199373,
                        63025520,
                        48129895,
                        51304566,
                        198022430,
                        191154276,
                        180915260,
                        171115067,
                        159138663,
                        146364022,
                        141213431,
                        16571,
                        155270560,
                        59373566)
export(lbingr, paste0(dir, "1kbin_time.bigWig"))

ggplot(ldfbin) +
	aes(x=start, y=total_usec) +
	geom_point() +
	facet_wrap(~ seqname)

ldfbin %>% arrange(desc(total_usec)) %>% View()

# BLACKLIST threshold
library(purple)
detach("package:purple", unload=TRUE)
library(purple)
library(RMySQL)
source("../gridss/libgridss.R")
library(tidyverse)
library(Biostrings)
db = dbConnect(MySQL(), dbname = "hmfpatients_sliced", host="localhost", port=3306, username="test", password="test")
callgr = query_structural_variants_for_sample_as_granges(db, query_structural_variants_samples(db))

time_portion = function(threshold) sum(as.numeric(lbingr[lbingr$total_usec > threshold]$total_usec)) / sum(as.numeric(lbingr$total_usec))
width_portion = function(threshold) sum(as.numeric(width(lbingr[lbingr$total_usec > threshold]))) / sum(as.numeric(width(lbingr)))
lost_calls = function(threshold) sum(overlapsAny(callgr, lbingr[lbingr$total_usec > threshold], ignore.strand=TRUE))

time_width_df = data.frame(threshold=10**((1:100)/10)) %>%
  mutate(time=Vectorize(time_portion)(threshold),
         width=Vectorize(width_portion)(threshold))
time_width_df = time_width_df %>%
  mutate(calls = Vectorize(lost_calls)(threshold))
ggplot(time_width_df) +
  aes(x=time * 100, y=width * 100, label=log10(threshold)) +
  geom_label() +
  coord_cartesian(xlim=c(0, 25), y=c(0, 1))

ggplot(time_width_df) +
  aes(x=time * 100, y=width * 100, label=log10(threshold)) +
  geom_label() +
  coord_cartesian(xlim=c(0, 25), y=c(0, 1))

ggplot(time_width_df) +
  aes(x=time * 100, y=calls, label=log10(threshold)) +
  geom_label() +
  coord_cartesian(xlim=c(0, 25), y=c(0, 500))


#dbNew = dbConnect(MySQL(), dbname = "hmfpatients_pilot")

theme_set(theme_bw())
lbingr66 <- lbingr[log10(lbingr$total_usec) > 6.6]

