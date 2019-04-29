#library(BiocManager)
#install("VariantAnnotation")
#install("dplyr")
#install("ggplot2")

library(dplyr)
library(GenomicRanges)
library(VariantAnnotation)
library(cowplot)
theme_set(theme_bw())

vcf_data_frame<- function(vcf) {
  vcf.rowRanges = rowRanges(vcf)
  vcf.info = info(vcf)
  vcf.alt = CharacterList(alt(vcf))
  vcf.alt[elementNROWS(vcf.alt) > 1L ] <- lapply(vcf.alt[elementNROWS(vcf.alt) > 1L ], paste0, collapse=",")
  
  vcf.df = data.frame(
    chr = seqnames(vcf), 
    pos = start(vcf), 
    ref = ref(vcf), 
    alt = as.character(vcf.alt),  
    filter = as.character(vcf.rowRanges$FILTER), 
    minorAllelePloidy = vcf.info$PURPLE_MAP, 
    ploidy = vcf.info$PURPLE_PLOIDY
    )
  
  vcf.df = vcf.df %>% mutate(alt = as.character(alt), type = ifelse(nchar(ref) == nchar(alt), "SNV", "INDEL"))
  
  return (vcf.df)
}

sample = "CPCT02050339"
sample = "CPCT02090044"

first.vcf = readVcf(file = paste0("~/hmf/analysis/deep/", sample, "T/", sample, "A.purple.somatic.vcf.gz"), 'hg19')
second.vcf = readVcf(file = paste0("~/hmf/analysis/deep/", sample, "T/", sample, "B.purple.somatic.vcf.gz"), 'hg19')
combined.vcf = readVcf(file = paste0("~/hmf/analysis/deep/", sample, "T/", sample, "C.purple.somatic.vcf.gz"), 'hg19')

first.df = vcf_data_frame(first.vcf)
second.df = vcf_data_frame(second.vcf)
combined.df = vcf_data_frame(combined.vcf)


first.df %>% filter(filter == 'PASS') %>% group_by(type) %>% dplyr::count()



scope <- function(df, vcf1, vcf2, vcf3) {
  ol = as.matrix(findOverlaps(vcf1, vcf2, type="any", select="all"))
  df$vcf2 <- F
  df[ol[, 1], "vcf2"] = T 
  df[ol[, 1], "vcf2Ploidy"] = info(vcf2)$PURPLE_PLOIDY[ol[, 2]] 
  
  ol = as.matrix(findOverlaps(vcf1, vcf3, type="any", select="all"))
  df$vcf3 <- F
  df[ol[, 1], "vcf3"] = T 
  df[ol[, 1], "vcf3Ploidy"] = info(vcf3)$PURPLE_PLOIDY[ol[, 2]] 
  
  df = df %>% mutate(
    scope = ifelse(vcf2, "vcf2", "Private"),
    scope = ifelse(vcf3, "vcf3", scope),
    scope = ifelse(vcf2&vcf3, "Shared", scope)) 
  
  return (df)
}

first.df = scope(first.df, first.vcf, second.vcf, combined.vcf) %>%
  mutate(
    scope = ifelse(scope == "vcf2", "Duplicate", scope), 
    scope = ifelse(scope == "vcf3", "Combined", scope))

second.df = scope(second.df, second.vcf, first.vcf, combined.vcf) %>%
  mutate(
    scope = ifelse(scope == "vcf2", "Production", scope), 
    scope = ifelse(scope == "vcf3", "Combined", scope))

combined.df = scope(combined.df, combined.vcf, first.vcf, second.vcf) %>%
  mutate(
    scope = ifelse(scope == "vcf2", "Production", scope), 
    scope = ifelse(scope == "vcf3", "Duplicate", scope))





ploidy_plot <- function (df) {
  ggplot(df %>% filter(filter == 'PASS'), aes(x = ploidy)) +
    geom_histogram(aes(fill = scope), alpha = 1,  binwidth = 0.1, color = "black",  position = "stack") + 
    scale_x_continuous(breaks = c(0:10), limits = c(-0.1, 3)) +
    theme(panel.grid.minor = element_blank(), axis.ticks = element_blank(), legend.title = element_blank()) +
    xlab("Ploidy") + ylab("Count") + ggtitle("Ploidy PDF")
}


p1 = ploidy_plot(first.df) + ggtitle(paste0(sample, "T Production"))
p2 = ploidy_plot(second.df) + ggtitle(paste0(sample, "T Duplicate"))
p3 = ploidy_plot(combined.df) + ggtitle(paste0(sample, "T Combined"))

p4 = ggplot(combined.df %>% filter(filter == 'PASS', !is.na(vcf2Ploidy), !chr %in% c('X', 'Y')), aes(x = ploidy, y = vcf2Ploidy)) +
  geom_point(alpha = 1, color = "black") +
  theme(panel.grid.minor = element_blank(), axis.ticks = element_blank(), legend.title = element_blank()) +
  xlab("Combined Ploidy ") + ylab("Production Ploidy")
p4


p5 = ggplot(combined.df %>% filter(filter == 'PASS', !is.na(vcf3Ploidy), !chr %in% c('X', 'Y')), aes(x = ploidy, y = vcf3Ploidy)) +
  geom_point(alpha = 1, color = "black") +
  theme(panel.grid.minor = element_blank(), axis.ticks = element_blank(), legend.title = element_blank()) +
  xlab("Combined Ploidy") + ylab("Duplicate Ploidy")

pTop = plot_grid(p1, p2, p3, nrow = 1)
pBottom = plot_grid(p4, p5,  nrow = 1)
pFinal = plot_grid(pTop, pBottom, nrow = 2)
pFinal
save_plot(paste0("/Users/jon/hmf/analysis/deep/",sample, "T.png"), pFinal, base_width = 17, base_height = 10)


combined.df %>% filter(filter == 'PASS') %>% group_by(scope) %>% dplyr::count()



##### SECOND ATTEMPT

summary <- function(sample) {
  first.vcf = readVcf(file = paste0("~/hmf/analysis/deep/", sample, "T/", sample, "A.purple.somatic.vcf.gz"), 'hg19')
  second.vcf = readVcf(file = paste0("~/hmf/analysis/deep/", sample, "T/", sample, "B.purple.somatic.vcf.gz"), 'hg19')
  combined.vcf = readVcf(file = paste0("~/hmf/analysis/deep/", sample, "T/", sample, "C.purple.somatic.vcf.gz"), 'hg19')
  
  first.df = vcf_data_frame(first.vcf)
  second.df = vcf_data_frame(second.vcf)
  combined.df = vcf_data_frame(combined.vcf)
  
  first.summary = first.df %>% filter(filter == 'PASS') %>% group_by(type) %>% dplyr::count() %>% mutate(source = "A", sample = sample)
  second.summary = first.df %>% filter(filter == 'PASS') %>% group_by(type) %>% dplyr::count() %>% mutate(source = "B", sample = sample)
  third.summary = first.df %>% filter(filter == 'PASS') %>% group_by(type) %>% dplyr::count() %>% mutate(source = "C", sample = sample)
  
  return (bind_rows(bind_rows(first.summary, second.summary), third.summary))
}


sample1 = summary('CPCT02050339') %>% mutate(sample = "Sample1")
sample2 = summary('CPCT02090044') %>% mutate(sample = "Sample2")
data = bind_rows(sample1, sample2)



plot_changes <- function(data, multiplier, percentageOffset, title) {
  colnames(data) <- c("sampleId", "prod", "down")
  totalAverage = data %>% ungroup() %>% summarise(prod = sum(prod), down = sum(down)) %>% mutate(change = (prod - down) / prod) %>% pull(change)
  tidyData = data %>% mutate(change = (prod - down) / prod) %>% gather(variable, value, 2, 3) %>% ungroup() %>% mutate(variable = ifelse(variable == "prod", "Normal Depth","Downsampled"))
  tidyData$average = totalAverage
  
  sampleFactor = tidyData %>% filter(variable == "Normal Depth") %>% arrange(-value) %>% pull(sampleId)
  variableFactor = c("Normal Depth","Downsampled")
  variableFactorColours = setNames(c("#6baed6", "#d94701"), variableFactor)
  tidyData = tidyData %>% mutate(sampleId = factor(sampleId, sampleFactor), variable = factor(variable, variableFactor))
  
  ggplot(data, aes(x = sampleId)) + 
    geom_bar(alpha = 0.8, aes(x = sampleId, y = value, fill = variable), stat = "identity", position = "dodge") +
    geom_point(aes(y = multiplier * (percentageOffset + change) )) +
    #geom_line(aes(y = multiplier * (percentageOffset + change), group = 1)) +
    geom_line(aes(y = multiplier * (percentageOffset + average)), group = 1, linetype = "dashed") + 
    xlab("Samples") + ylab("Count") + ggtitle(title) + 
    scale_y_continuous(expand = c(0,0), sec.axis = sec_axis(~./multiplier - percentageOffset, name = " % difference", labels = percent)) +
    scale_fill_manual(name = "", values = variableFactorColours) +
    theme(panel.border = element_blank(), axis.ticks = element_blank(), axis.text.x = element_blank(), legend.position = c(0.7,0.9), panel.grid.major.x = element_blank())
}
