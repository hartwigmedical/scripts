get_sig_data<-function(signatures, contribution, sigNames, sampleNames) {

  # get contributions by sample in percentage terms (ie each sample split across the signatures)
  contributionPercents = apply(contribution, 2, function(x) x/sum(x))
  transposedPerc = t(contributionPercents) %>% as.data.frame()
  colnames(transposedPerc) <- sigNames
  sigPercentsBySample = cbind(data.frame(SampleId = sampleNames), transposedPerc)
  rownames(sigPercentsBySample) <- NULL

  # get contributions by sample in raw count of SVs
  transposedContrib = t(contribution) %>% as.data.frame()
  colnames(transposedContrib) <- sigNames

  # use signature counts to turn contributions into counts
  sigCountsBySample = transposedContrib
  for (i in 1:ncol(sigCountsBySample))
  {
    sigCountsBySample[i] = round(sigCountsBySample[i]* colSums(signatures)[i],0)
  }

  sigCountsBySample = cbind(data.frame(SampleId = sampleNames), sigCountsBySample)
  rownames(sigCountsBySample) <- NULL

  # now put percentages and raw counts together, and split columns into rows to help with sig-agnostic aggregation
  sampleSigPercData = sigPercentsBySample %>% gather(SigName, SigPercent, -SampleId) %>% arrange(SampleId, SigName)
  sampleSigPercData$PercBucket = round(sampleSigPercData$SigPercent/0.1)*0.1

  # do same for counts then add together
  sampleSigCountData = sigCountsBySample %>% gather(SigName, SvCount, -SampleId) %>% arrange(SampleId, SigName)

  # ensure both are order by sample and sig
  sampleSigData = cbind(sampleSigPercData, SvCount = sampleSigCountData$SvCount)

  return (sampleSigData)
}

get_sig_stats<-function(sampleSigData) {

  # key stats per signature
  sigStats = (sampleSigData %>% group_by(SigName)
              %>% summarise(SampleCount=sum(SvCount>0),
                            SamplePerc=round(sum(SvCount>0)/n_distinct(SampleId),2),
                            SvCount=sum(SvCount),
                            MaxPercent=round(max(SigPercent),3),
                            AvgPercent=round(sum(ifelse(SigPercent>0.001,SigPercent,0))/sum(SigPercent>0),3),
                            PercGT75=sum(SigPercent>0.75),
                            Perc50_75=sum(SigPercent>0.5&SigPercent<=0.75),
                            Perc25_50=sum(SigPercent>0.25&SigPercent<=0.5),
                            Perc10_25=sum(SigPercent>0.1&SigPercent<=0.25),
                            Perc5_10=sum(SigPercent>0.05&SigPercent<=0.1),
                            PercLT5=sum(SigPercent>0.001&SigPercent<=0.05))
              %>% arrange(-SampleCount))

  return (sigStats)
}


plot_sig_samples<-function(sampleSigData, cancerType)
{
  cancerSigData = sampleSigData

  if(cancerType != "") {
    cancerSigData = cancerSigData %>% filter(CancerType==cancerType)
  }

  sigColours = c("#ff994b","#463ec0","#d10073","#996ffb","#68b1c0","#e34bd9","#106b00","#8a392f","#98d76a","#6b3a9d","#d5c94e","#c6c0fb",
                 "#ff862c","#31528d","#d7003a","#323233","#ff4791","#01837a", "#ff748a","#777700","#ff86be")

  cancerSampleSigData = cancerSigData %>% arrange(-SampleSvCount, SampleId) %>% select('SampleId', 'SigName', 'SvCount')

  sigCancerPlots = list()
  plotIndex = 1

  cancerSigStats = svnmf::get_sig_stats(cancerSigData)

  if(nrow(cancerSigStats) > 0) {

    if(cancerType == "")
    {
      title = "Signature SV Counts"
    }
    else
    {
      title = paste("Signature SV Counts for ", cancerType, sep="")
    }

    sigStatsPlot = (ggplot(data = cancerSigStats, aes(x = SigName, y = SvCount, group = 1), fill = SigName)
                    + geom_bar(stat = "identity", colour = "black", size = 0.2)
                    + theme(axis.text.x = element_text(angle = 90, hjust = 1))
                    + ylab("SV Count") + xlab("Signature") + ggtitle(title)
                    + theme(legend.position="none"))

    axisRatio = max(cancerSigStats$SvCount) / max(cancerSigStats$SampleCount) * 0.6

    sigStatsPlot <- (sigStatsPlot + geom_line(aes(y = SampleCount*axisRatio, color = "red"))
                     + scale_y_continuous(sec.axis = sec_axis(~.*(1/axisRatio), name = "Sample Count")))

    sigCancerPlots[[plotIndex]] <- sigStatsPlot
    plotIndex = plotIndex + 1
  }

  if(nrow(cancerSampleSigData) > 0) {

    # only plot 50 samples at a time
    numSamples = n_distinct(cancerSampleSigData$SampleId)
    numSigs = n_distinct(cancerSampleSigData$SigName)
    samplesPerPlot = 50
    rowsPerPlot = samplesPerPlot * numSigs
    plotCount = ceiling(numSamples/samplesPerPlot)

    for (n in 1:plotCount) {
      rowStart = ((n-1) * rowsPerPlot + 1)
      rowEnd = min((n * rowsPerPlot), numSamples * numSigs)

      if(cancerType == "")
      {
        title = "Sig SV Counts by Sample"
      }
      else
      {
        title = paste("Sig SV Counts by Sample for ", cancerType, sep="")
      }

      sampleSigPlot <- (ggplot(cancerSampleSigData[rowStart:rowEnd,], aes(x = reorder(SampleId, -SvCount), y = SvCount, fill = SigName))
                        + geom_bar(stat = "identity", colour = "black")
                        + labs(x = "", y = "SV Count by Sample")
                        + scale_fill_manual(values = sigColours)
                        + theme_bw() + theme(panel.grid.minor.x = element_blank(), panel.grid.major.x = element_blank())
                        + theme(panel.grid.minor.y = element_blank(), panel.grid.major.y = element_blank())
                        + theme(axis.text.x = element_text(angle = 90, hjust = 1,size=7)))

      if(n == 1)
      {
        sampleSigPlot <- sampleSigPlot + ggtitle(title)
      }
      if(n > 1)
      {
        # remove legend after the first plot
        sampleSigPlot <- sampleSigPlot + theme(legend.position="none")
      }

      sigCancerPlots[[plotIndex]] <- sampleSigPlot

      if(plotIndex >=4)
      {
        multiplot(plotlist = sigCancerPlots, cols = 2)
        sigCancerPlots = list()
        plotIndex = 1
      }
      else
      {
        plotIndex = plotIndex + 1
      }
    }

    if(plotIndex > 1)
    {
      # now print all plots for this cancer type
      multiplot(plotlist = sigCancerPlots, cols = 2)
    }
  }
}

plot_top_n_samples_by_sig<-function(sampleSigData, sigNames, topN = 50)
{

  sigColours = c("#ff994b","#463ec0","#d10073","#996ffb","#68b1c0","#e34bd9","#106b00","#8a392f","#98d76a","#6b3a9d","#d5c94e","#c6c0fb",
                 "#ff862c","#31528d","#d7003a","#323233","#ff4791","#01837a", "#ff748a","#777700","#ff86be")

  sigSamplePlots = list()
  plotIndex = 1

  # merge cancer type with sampleId
  sampleSigData = unite(sampleSigData, "SampleId", SampleId, CancerType, sep='_')

  for(sigName in sigNames)
  {
    topNSamplesBySig = head(sampleSigData %>% filter(SigName==sigName) %>% arrange(-SvCount),topN)
    # View(topNSamplesBySig)

    # now grab all sig data for these top-N samples
    topNSampleSigData = sampleSigData %>% filter(SampleId %in% topNSamplesBySig$SampleId)
    # View(top30SampleSigData)

    title = paste("Top Samples for Signature ", sigName, sep="")

    sampleSigPlot <- (ggplot(topNSampleSigData, aes(x = reorder(SampleId, -SvCount), y = SvCount, fill = SigName))
                      + geom_bar(stat = "identity", colour = "black")
                      + labs(x = "", y = "SV Count by Sample")
                      + scale_fill_manual(values = sigColours)
                      + theme_bw() + theme(panel.grid.minor.x = element_blank(), panel.grid.major.x = element_blank())
                      + theme(panel.grid.minor.y = element_blank(), panel.grid.major.y = element_blank())
                      + theme(axis.text.x = element_text(angle = 90, hjust = 1,size=7))
                      + ggtitle(title))

    if(plotIndex > 1)
    {
      # remove legend after the first plot
      sampleSigPlot <- sampleSigPlot + theme(legend.position="none")
    }

    sigSamplePlots[[plotIndex]] <- sampleSigPlot

    if(plotIndex >=4)
    {
      multiplot(plotlist = sigSamplePlots, cols = 2)
      sigSamplePlots = list()
      plotIndex = 1
    }
    else
    {
      plotIndex = plotIndex + 1
    }
  }

  if(plotIndex > 1)
  {
    # now print all plots for this cancer type
    multiplot(plotlist = sigSamplePlots, cols = 2)
  }
}
