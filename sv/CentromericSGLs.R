
dataDir = '~/data/sv/sgl_blat'


# Plot 1: Destination of centromeric SGLs
write.csv(sglDestChrs,'~/data/sv/sgl_blat/sgl_destination_chromosomes.csv',row.names = F,quote = F)
sglDestChrs = read.csv(paste(dataDir,'sgl_destination_chromosomes.csv',sep='/'))

print(ggplot(sglDestChrs, aes(x=reorder(Chromosome,ChrIndex), y=n))
      + geom_bar(stat = "identity", colour = "black")
      + labs(x='Chromosome', y='# Centromeric SGLs', title='Desintation for Centromeric SGLs'))


# Plot 2: Frequency of SVs ending on each chromosome and arm by whether a local or translocation SV
write.csv(sglLocalVsTrans,'~/data/sv/sgl_blat/sgl_local_vs_trans.csv',row.names = F,quote = F)
sglLocalVsTrans = read.csv(paste(dataDir,'sgl_local_vs_trans.csv',sep='/'))

print(ggplot(sglLocalVsTrans, aes(x=reorder(CentroChrText,-Total), y=Percent, fill=SvType))
      + geom_bar(stat = "identity", colour = "black")
      + labs(x='Centromere Chromosome (total SGLs)', y='% of Centromeric SGLs per Chromosome',title='Proportion of Local vs Translocation Centromeric SGLs by Chromosome')
      + theme_bw() + theme(panel.grid.minor.x = element_blank(), panel.grid.major.x = element_blank())
      + theme(panel.grid.minor.y = element_blank(), panel.grid.major.y = element_blank())
      + theme(axis.text.x = element_text(angle = 90, hjust = 1,size=10)))



# Plot 3: Location of local SGLs by chromosome
write.csv(sglOrigPositions,'~/data/sv/sgl_blat/sgl_originating_positions.csv',row.names = F,quote = F)
sglOrigPositions = read.csv(paste(dataDir,'sgl_originating_positions.csv',sep='/'))

print(ggplot(sglOrigPositions, aes(x=reorder(SourcePosLabel,PosBucket), y=n))
      + geom_bar(stat = "identity", colour = "black")
      + facet_wrap(~Chromosome)
      + labs(x='Chromosomal Position', y='# SGLs', title = "Originating position of Local Centromeric SGLs by Chromosome")
      + theme(axis.text.x=element_blank(),axis.ticks.x=element_blank()))


# Plot 4: Centromeric SGLs vs Arms with CN Gain
write.csv(sglFbCentroSummary,'~/data/sv/sgl_blat/sgl_arm_cn_gain_summary.csv',row.names = F,quote = F)
sglFbCentroSummary = read.csv(paste(dataDir,'sgl_arm_cn_gain_summary.csv',sep='/'))

print(ggplot(sglFbCentroSummary %>% filter(CentroGain>0|SglCentro>0|FBArms>0), aes(x=CentroGain,y=SglCentro))
      + geom_point(position="jitter")
      + geom_smooth(,method=lm,se=FALSE, fullrange=F)
      + labs(x='Arms with Centromeric Copy Number Change', y='Arms with Centromeric SGLs',title = "Chromosomal Arms per Sample with Centromeric Copy Number Change vs Centromeric SGLs"))

# Plot 5: Rate of CN Gain vs Centromeric SGLs
sglRateOfGain = read.csv(paste(dataDir,'sgl_rate_of_gain.csv',sep='/'))

print(ggplot(sglRateOfGain, aes(x=reorder(Chromosome,ChrIndex),y=Rate,fill=RateOfCentromericGain))
      + geom_bar(position="dodge",stat='identity')
      + labs(x='Centromere Chromosome', y='% of Samples with Centromeric Gain',title='Rate of Centromeric CN Gain with and without Centromeric SGLs'))



