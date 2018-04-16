library(tidyverse)

fdf = read.csv("D://hartwig//CPCT02100145.gridss.assembly.bam.events.csv", header=FALSE)
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




# binned timing
ldfbin <- ldf %>% 
	mutate(bin=round(start, -4)) %>%
	group_by(chunk, direction, seqname, bin) %>%
	summarise(
		total_usec=sum(usec),
		max_records=max(records))

ggplot(ldfbin) +
	aes(x=bin, y=total_usec) +
	geom_point() + 
	facet_wrap(~ seqname)

