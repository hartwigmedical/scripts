get_sig_colours<-function(sigCount = 10)
{
  # original COSMIC
  # allColours = c("#ff994b", "#463ec0", "#88c928","#996ffb","#68b1c0","#e34bd9","#106b00","#d10073","#98d76a",
  #              "#6b3a9d","#d5c94e","#0072e2","#ff862c","#31528d","#d7003a","#323233","#ff4791","#01837a",
  #              "#ff748a","#777700","#ff86be","#4a5822","#ffabe4","#6a4e03","#c6c0fb","#ffb571","#873659",
  #              "#dea185","#a0729d","#8a392f")

  allColours = c("#ff994b", "#463ec0", "#88c928", "#996ffb", "#68b1c0", "#e34bd9", "#106b00", "#d10073", "#98d76a",
                 "#6b3a9d", "#d5c94e", "#0072e2", "red3", "#31528d", "aquamarine", "tomato2", "#ff4791", "#01837a",
                 "#ff748a", "#777700", "#ff86be", "#4a5822", "#ffabe4", "#6a4e03", "#c6c0fb", "#ffb571", "#873659",
                 "#ff494b", "#464ec0", "#885928", "#991ffb", "#6811c0", "#e38bd9", "#101b00", "#d12073", "#98176a",
                 "#6b4a9d", "#d5494e", "#0052e2", "#ff162c", "#31128d", "#d7803a", "#421233", "#ff2791", "#01137a",
                 "#ff448a", "#774700", "#ff56be", "#4a1822", "#ff1be4", "#6a8e03", "#c610fb", "#ff2571", "#871659",
                 "#ff294b", "#465ec0", "#886928", "#992ffb", "#6821c0", "#e37bd9", "#102b00", "#d13073", "#98276a",
                 "#6b2a9d", "#d5594e", "#0062e2", "#ff262c", "#31228d", "#47703a", "#422233", "#ff3791", "#01237a",
                 "#ff248a", "#775700", "#ff66be", "#4a2822", "#ff2be4", "#6a7e03", "#c620fb", "#ff3571", "#872659",
                 "#dea185", "#a0629d", "#8a792f", "#4a3822", "#ff3be4", "#6a6e03", "#c630fb", "#ff4571", "#877659")

  sigColours = c()

  for(i in 1:sigCount)
  {
    sigColours[i] = allColours[i]
  }

  return (sigColours)
}

get_unique_colours<-function()
{
  colourSet = c("cornsilk3", "hotpink", "darkorange", "seagreen3", "tomato3", "thistle2", "steelblue2", "darkgreen", "indianred", "honeydew2",
                "turquoise3", "lightpink2", "goldenrod2", "darkslateblue", "yellowgreen", "wheat2", "violetred2", "ivory3", "coral1", "springgreen2")

  return (colourSet)
}

get_base_colours<-function()
{
  baseColours = c("yellow", "blue", "green", "red", "orange", "purple", "pink", "brown", "darkgreen", "deepskyblue", "tan")
  # baseColours = c("yellow", "blue", "green", "red", "orange", "purple", "pink", "brown")
  return (baseColours)
}

get_base_colour_extensions<-function(baseColour, extensionCount)
{
  if(baseColour == "yellow")
    extnColours = c("lightgoldenrod1", "lightgoldenrod2", "lightgoldenrod3", "lightgoldenrod4")
  else if(baseColour == "blue")
    extnColours = c("royalblue1", "royalblue2", "royalblue3", "royalblue4")
  else if(baseColour == "green")
    extnColours = c("palegreen1", "palegreen2", "palegreen3", "palegreen4")
  else if(baseColour == "red")
    extnColours = c("firebrick1", "firebrick2", "firebrick3", "firebrick4")
  else if(baseColour == "orange")
    extnColours = c("orangered1", "orangered2", "orangered3", "orangered4")
  else if(baseColour == "purple")
    extnColours = c("magenta1", "magenta2", "magenta3", "magenta4")
  else if(baseColour == "pink")
    extnColours = c("deeppink1", "deeppink2", "deeppink3", "deeppink4")
  else if(baseColour == "brown")
    extnColours = c("chocolate1", "chocolate2", "chocolate3", "chocolate4")
  else if(baseColour == "darkgreen")
    extnColours = c("darkolivegreen1", "darkolivegreen2", "darkolivegreen3", "darkolivegreen4")
  else if(baseColour == "deepskyblue")
    extnColours = c("cyan1", "cyan2", "cyan3", "cyan4")
  else if(baseColour == "tan")
    extnColours = c("tan1", "tan2", "tan3", "tan4")

  return (extnColours)
}

strip_multi_bg_colours<-function(sigColours, bgSigCount)
{
  # the background colour is repeated X times, strip back so it is only in there once
  newSigColours = c()
  sigCount = length(sigColours)

  # newSigColours[1] = sigColours[1]

  for(i in bgSigCount:sigCount)
  {
    # so if say BG count = 20, the 1st colour will be set to the 20th
    newSigColours[i-bgSigCount+1] = sigColours[i]
  }

  return (newSigColours)
}

get_signame_list<-function(sigCount, asStrings, zeroBased=F)
{
  sigs = c()
  for(i in 1:sigCount)
  {
    index = ifelse(zeroBased,i-1,i)

    if(asStrings & i <= 9)
    {
      sigs[i] = paste("0", index, sep='')
    }
    else
    {
      sigs[i] = index
    }
  }

  return (sigs)
}

# library(pracma)

add_missing_buckets<-function(matrixData, definedBuckets)
{
  sampleCount = ncol(matrixData)-1
  missingBuckets = definedBuckets %>% filter(!Bucket %in% matrixData$Bucket)

  for(missingBucket in missingBuckets$Bucket)
  {
    print(paste("adding missing bucket: ", missingBucket, sep=''))

    rowIndex = nrow(matrixData)+1
    matrixData[rowIndex,1] = missingBucket

    for(i in 1:sampleCount)
    {
      matrixData[rowIndex,i+1] = 0
    }
  }

  return (matrixData)
}

apply_signatures<-function(matrixData, signatures)
{
  n_samples = dim(matrixData)[2]
  n_signatures = dim(signatures)[2]
  lsq_contribution = matrix(NA, nrow = n_signatures, ncol = n_samples)

  for(i in 1:ncol(matrixData))
  {
    # print(paste("calc LSQ for i=", i, sep=''))
    y = matrixData[, i, drop=F]
    y = apply(y, 1, function(x) as.numeric(x) * 1.0) # to force it to a numeric type
    lsq = calc_lsqnonneg(signatures, y)
    lsq_contribution[, i] = lsq$x
    # lsq_reconstructed[, i] = signatures %*% as.matrix(lsq$x)
  }
  sample_names = colnames(matrixData)
  signature_names = colnames(signatures)
  mut_type_names = rownames(signatures)
  colnames(lsq_contribution) = sample_names
  rownames(lsq_contribution) = signature_names
  #colnames(lsq_reconstructed) = sample_names
  #rownames(lsq_reconstructed) = mut_type_names
  # res = list(lsq_contribution, lsq_reconstructed)
  # names(res) = c("contribution", "reconstructed")

  return(lsq_contribution)
}

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
    sigCountsBySample[i] = sigCountsBySample[i]* colSums(signatures)[i]
  }

  sigCountsBySample = cbind(data.frame(SampleId = sampleNames), sigCountsBySample)
  rownames(sigCountsBySample) <- NULL

  # now put percentages and raw counts together, and split columns into rows to help with sig-agnostic aggregation
  sampleSigPercData = sigPercentsBySample %>% gather(SigName, SigPercent, -SampleId) %>% arrange(SampleId, SigName)
  sampleSigPercData$PercBucket = round(sampleSigPercData$SigPercent/0.1)*0.1

  # do same for counts then add together
  sampleSigCountData = sigCountsBySample %>% gather(SigName, Count, -SampleId) %>% arrange(SampleId, SigName)

  # ensure both are order by sample and sig
  sampleSigData = cbind(sampleSigPercData, Count = sampleSigCountData$Count)

  return (sampleSigData)
}

append_residuals<-function(contribution, signatures, matrixData, bucketNames, sampleSigData)
{
  residuals = calc_contrib_sample_residuals(contribution, signatures, matrixData, bucketNames)

  sampleResiduals = residuals %>% select(SampleId,UnallocCount)
  colnames(sampleResiduals) <- c("SampleId","Count")

  sampleResiduals$SigName = "Unalloc"
  sampleResiduals$SigPercent = 0
  sampleResiduals$PercBucket = 0
  sampleSigData2 = rbind(sampleSigData, sampleResiduals %>% select(SampleId,SigName,SigPercent,PercBucket,Count))

  sampleResiduals = residuals %>% select(SampleId,ExcessCount)
  colnames(sampleResiduals) <- c("SampleId","Count")

  sampleResiduals$SigName = "Excess"
  sampleResiduals$SigPercent = 0
  sampleResiduals$PercBucket = 0
  sampleResiduals$Count = sampleResiduals$Count * -1.0
  sampleSigData2 = rbind(sampleSigData2, sampleResiduals %>% select(SampleId,SigName,SigPercent,PercBucket,Count))

  return (sampleSigData2)
}

calc_sig_alloc_sample_residuals<-function(sigSampleBucketCounts, matrixData, bucketNamesIndexed)
{
  sampleSigBucketTotals = sigSampleBucketCounts %>% group_by(SampleId, Bucket) %>% summarise(SigCount=sum(Count))

  # calculate residuals by lining up the actual vs the sig counts per sample per bucket
  sampleBucketCounts = matrixData
  sampleBucketCounts = cbind(bucketNamesIndexed$BucketIndex, sampleBucketCounts)
  gatherIndex = ncol(sampleBucketCounts)
  sampleBucketCounts2 = gather(as.data.frame(sampleBucketCounts), "SampleId", "ActualCount", 2:gatherIndex)
  colnames(sampleBucketCounts2) = c("Bucket", "SampleId", "ActualCount")

  sampleBucketFitVsActuals = merge(sampleBucketCounts2, sampleSigBucketTotals, by=c("SampleId","Bucket"), all=T)

  sampleBucketFitVsActuals$ResidualExcess = ifelse(sampleBucketFitVsActuals$SigCount > sampleBucketFitVsActuals$ActualCount, sampleBucketFitVsActuals$SigCount - sampleBucketFitVsActuals$ActualCount, 0)
  sampleBucketFitVsActuals$ResidualUnalloc = ifelse(sampleBucketFitVsActuals$ActualCount > sampleBucketFitVsActuals$SigCount, sampleBucketFitVsActuals$ActualCount - sampleBucketFitVsActuals$SigCount, 0)
  sampleBucketFitVsActuals$ResidualDiff = abs(sampleBucketFitVsActuals$ActualCount-sampleBucketFitVsActuals$SigCount)

  sampleResidualData = (sampleBucketFitVsActuals %>% group_by(SampleId)
                        %>% summarise(ActualCount=sum(ActualCount),
                                      SigCount=sum(SigCount),
                                      ResidualTotal=round(sum(ResidualDiff),2),
                                      ExcessCount=round(sum(ResidualExcess),2),
                                      UnallocCount=round(sum(ResidualUnalloc),2),
                                      ResidualPerc=round(sum(ResidualDiff)/sum(ActualCount),2)))

  sampleResidualData$ResidualTotalManual = sampleResidualData$ActualCount - sampleResidualData$SigCount

  return (sampleResidualData)
}

calc_contrib_sample_residuals<-function(contributions, signatures, matrixData, bucketNames)
{
  sampleBucketFit = signatures %*% contributions

  sampleBucketFit = cbind(bucketNames,sampleBucketFit)
  sampleBucketCounts = cbind(bucketNames,matrixData)

  gatherIndex = ncol(sampleBucketCounts)
  sampleBucketCounts2 = gather(sampleBucketCounts, "SampleId", "ActualCount", 2:gatherIndex)
  sampleBucketFit2 = gather(sampleBucketFit, "SampleId", "FitCount", 2:gatherIndex)

  sampleBucketFitVsActuals = merge(sampleBucketCounts2, sampleBucketFit2, by=c("SampleId","Bucket"), all=T)

  sampleBucketFitVsActuals$ResidualExcess = ifelse(sampleBucketFitVsActuals$FitCount > sampleBucketFitVsActuals$ActualCount, sampleBucketFitVsActuals$FitCount - sampleBucketFitVsActuals$ActualCount, 0)
  sampleBucketFitVsActuals$ResidualUnalloc = ifelse(sampleBucketFitVsActuals$ActualCount > sampleBucketFitVsActuals$FitCount, sampleBucketFitVsActuals$ActualCount - sampleBucketFitVsActuals$FitCount, 0)
  sampleBucketFitVsActuals$ResidualDiff = abs(sampleBucketFitVsActuals$ActualCount-sampleBucketFitVsActuals$FitCount)

  sampleResidualData = (sampleBucketFitVsActuals %>% group_by(SampleId)
                        %>% summarise(ActualCount=sum(ActualCount),
                                      FitCount=sum(FitCount),
                                      ResidualTotal=round(sum(ResidualDiff),2),
                                      ExcessCount=round(sum(ResidualExcess),2),
                                      UnallocCount=round(sum(ResidualUnalloc),2),
                                      ResidualPerc=round(sum(ResidualDiff)/sum(ActualCount),2)))

  sampleResidualData$ResidualTotalManual = sampleResidualData$ActualCount - sampleResidualData$FitCount

  return (sampleResidualData)
}


get_sig_stats<-function(sampleSigData) {

  # key stats per signature
  sigStats = (sampleSigData %>% group_by(SigName)
              %>% summarise(SampleCount=sum(Count>0),
                            SamplePerc=round(sum(Count>0)/n_distinct(SampleId),2),
                            Count=round(sum(Count),0),
                            MaxPercent=round(max(SigPercent),3),
                            AvgPercent=round(sum(ifelse(SigPercent>0.001,SigPercent,0))/sum(SigPercent>0),3),
                            PercGT75=sum(SigPercent>0.75),
                            Perc50_75=sum(SigPercent>0.5&SigPercent<=0.75),
                            Perc25_50=sum(SigPercent>0.25&SigPercent<=0.5),
                            Perc10_25=sum(SigPercent>0.1&SigPercent<=0.25),
                            Perc5_10=sum(SigPercent>0.05&SigPercent<=0.1),
                            PercLT5=sum(SigPercent>0.001&SigPercent<=0.05))
              %>% arrange(SigName))

  return (sigStats)
}


plot_cancer_sigs<-function(sampleSigData, sigColours, varType = "SV")
{
  cancerSigData = sampleSigData %>% group_by(CancerType,SigName) %>% summarise(Count=sum(Count))
  cancerTotals = sampleSigData %>% group_by(CancerType) %>% summarise(TotalCount=sum(Count))
  cancerSigData = merge(cancerSigData, cancerTotals, by.x="CancerType", by.y="CancerType", all.x=TRUE)
  cancerSigData$Percent = cancerSigData$Count / cancerSigData$TotalCount

  cancerSigPlot <- (ggplot(cancerSigData, aes(x = CancerType, y = Percent, fill = SigName))
                    + geom_bar(stat = "identity", colour = "black")
                    + ggtitle("Signatures % by Cancer")
                    + labs(x = "", y = paste(varType, " % by Signature", sep=''))
                    + scale_fill_manual(values = sigColours)
                    + theme_bw() + theme(panel.grid.minor.x = element_blank(), panel.grid.major.x = element_blank())
                    + theme(panel.grid.minor.y = element_blank(), panel.grid.major.y = element_blank())
                    + theme(axis.text.x = element_text(angle = 90, hjust = 1,size=15)))

  print(cancerSigPlot)
}

get_sig_plot_legend<-function(sampleSigData, sigColours)
{
  sampleSigPlot <- (ggplot(sampleSigData, aes(x = SampleId, y = Count, fill = SigName))
                    + geom_bar(stat = "identity", colour = "black")
                    + scale_fill_manual(values = sigColours)
                    + theme_bw() + theme(panel.grid.minor.x = element_blank(), panel.grid.major.x = element_blank()))

  legend = get_legend(sampleSigPlot)
  return (legend)
}

plot_sig_samples<-function(sampleSigData, cancerType, sigColours, varType = "SV", maxPlots = 0, dropLegend = F)
{
  cancerSigData = sampleSigData

  if(cancerType != "")
  {
    cancerSigData = cancerSigData %>% filter(CancerType==cancerType)

    # filter out irrelevant background sig data if present
    countsByBGSigType = cancerSigData %>% filter(grepl("BG_",SigName)) %>% group_by(SigName) %>% summarise(Count=sum(Count)) %>% arrange(-Count)

    if(nrow(countsByBGSigType) > 0 && head(countsByBGSigType,1)$Count > 0)
    {
      topBGType = head(countsByBGSigType,1)$SigName
      cancerSigData = cancerSigData %>% filter((!grepl("BG_",SigName)|(SigName==topBGType)))
    }
  }

  cancerSampleSigData = cancerSigData %>% arrange(-SampleCount, SampleId) %>% select('SampleId','SigName','Count','SampleCount')

  sigCancerPlots = list()
  plotIndex = 1

  cancerSigStats = get_sig_stats(cancerSigData)

  if(nrow(cancerSigStats) > 0) {

    if(cancerType == "")
    {
      title = "Signature Counts All Cancer Types"
    }
    else
    {
      title = paste("Signature Counts for ", cancerType, sep="")
    }

    sigStatsPlot = (ggplot(data = cancerSigStats, aes(x = SigName, y = Count, group = 1), fill = SigName)
                    + geom_bar(stat = "identity", colour = "black", size = 0.2)
                    + theme(axis.text.x = element_text(angle = 90, hjust = 1))
                    + ylab(paste(varType, " Count", sep='')) + xlab("Signature") + ggtitle(title)
                    + theme(legend.position="none"))

    axisRatio = max(cancerSigStats$Count) / max(cancerSigStats$SampleCount) * 0.6

    sigStatsPlot <- (sigStatsPlot + geom_line(aes(y = SampleCount*axisRatio, color = "red"))
                     + scale_y_continuous(sec.axis = sec_axis(~.*(1/axisRatio), name = "Sample Count")))

    sigCancerPlots[[plotIndex]] <- sigStatsPlot
    plotIndex = plotIndex + 1
  }

  if(nrow(cancerSampleSigData) > 0)
  {
    # log the top N samples separately if they significantly higher counts
    cancerSampleData = cancerSigData %>% group_by(SampleId) %>% summarise(SampleCount=first(SampleCount)) %>% arrange(-SampleCount)

    topNIndex = 6

    logTopNSamples = F
    if(nrow(cancerSampleData) > topNIndex)
    {
      maxCount = max(cancerSampleData$SampleCount)
      nthCount = nth(cancerSampleData$SampleCount, topNIndex)

      if(maxCount > 3 * nthCount)
      {
        logTopNSamples = T
      }
    }

    # only plot 50 samples at a time
    numSamples = nrow(cancerSampleData)
    numSigs = n_distinct(cancerSampleSigData$SigName)
    samplesPerPlot = 50
    rowsPerPlot = samplesPerPlot * numSigs
    # plotCount = ceiling(numSamples/samplesPerPlot)

    rowEnd = 0
    sampleEnd = 0
    maxRows = nrow(cancerSampleSigData)
    plotCount = 0

    print(paste("cancer=", cancerType, ", numSamples=", numSamples, ", maxRows=", maxRows, ", logTopN=", logTopNSamples, sep=''))

    # while(rowEnd < maxRows)
    while(sampleEnd < numSamples)
    {
      if(plotCount == 0 & logTopNSamples)
      {
        sampleStart = 1
        sampleEnd = sampleStart + topNIndex - 1
        rowStart = 1
        rowEnd = min(rowStart + topNIndex*numSigs - 1, maxRows)
      }
      else
      {
        sampleStart = sampleEnd + 1
        sampleEnd = sampleStart + samplesPerPlot - 1
        rowStart = rowEnd + 1
        rowEnd = min(rowStart + rowsPerPlot - 1, maxRows)

        # combine stragglers onto the last graph
        if(rowEnd < maxRows & maxRows - rowEnd <= rowsPerPlot * 0.4)
        {
          # print(paste(plotCount, ": cancer=", cancerType, ", rowEnd=", rowEnd, " close to maxRows=", maxRows, sep=''))
          rowEnd = maxRows
        }

        if(sampleEnd < numSamples & numSamples - sampleEnd <= samplesPerPlot * 0.4)
        {
          sampleEnd = numSamples
        }
      }

      plotSampleSet = cancerSampleData[sampleStart:sampleEnd,]

      print(paste("plotCount= ", plotCount, ", sampleStart=", sampleStart, ", sampleEnd=", sampleEnd, sep=''))
      # print(paste("plotCount= ", plotCount, ", rowStart=", rowStart, ", rowEnd=", rowEnd, sep=''))

      if(cancerType == "")
      {
        title = "Sig Counts by Sample"
      }
      else
      {
        title = paste("Sig Counts by Sample for ", cancerType, sep="")
      }

      plotDataSet = cancerSampleSigData %>% filter(SampleId %in% plotSampleSet$SampleId)
      sampleSigPlot <- (ggplot(plotDataSet, aes(x = reorder(SampleId, -SampleCount), y = Count, fill = SigName))
                        + geom_bar(stat = "identity", colour = "black")
                        + labs(x = "", y = paste(varType, " Count by Sample", sep=''))
                        + scale_fill_manual(values = sigColours)
                        + theme_bw() + theme(panel.grid.minor.x = element_blank(), panel.grid.major.x = element_blank())
                        + theme(panel.grid.minor.y = element_blank(), panel.grid.major.y = element_blank())
                        + theme(axis.text.x = element_text(angle = 90, hjust = 1,size=7)))

      if(plotCount == 0)
      {
        # sampleSigPlot <- sampleSigPlot + ggtitle(title)

        if(dropLegend)
          sampleSigPlot <- sampleSigPlot + theme(legend.position="none")
      }
      else
      {
        # remove legend after the first plot
        sampleSigPlot <- sampleSigPlot + theme(legend.position="none")
      }

      sigCancerPlots[[plotIndex]] <- sampleSigPlot
      plotCount = plotCount + 1

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

      if(maxPlots > 0 & plotIndex > maxPlots)
        break
    }

    if(plotIndex > 1)
    {
      # now print all plots for this cancer type
      multiplot(plotlist = sigCancerPlots, cols = 2)
    }
  }
}

plot_top_n_samples_by_sig<-function(sampleSigData, sigNames, topN = 50, sigColours, varType = "SV", dropSigLegend = F)
{
  sigSamplePlots = list()
  plotIndex = 1

  # merge cancer type with sampleId
  sampleSigData = unite(sampleSigData, "SampleId", SampleId, CancerType, sep='_')

  for(sigName in sigNames)
  {
    topNSamplesBySig = head(sampleSigData %>% filter(SigName==sigName&Count>0) %>% arrange(-Count),topN)

    # now grab all sig data for these top-N samples
    topNSampleSigData = sampleSigData %>% filter(SampleId %in% topNSamplesBySig$SampleId)

    title = paste("Top Samples for Signature ", sigName, sep="")

    sampleSigPlot <- (ggplot(topNSampleSigData, aes(x = reorder(SampleId, -Count), y = Count, fill = SigName))
                      + geom_bar(stat = "identity", colour = "black")
                      + labs(x = "", y = paste(varType, " Count by Sample", sep=''))
                      + scale_fill_manual(values = sigColours)
                      + theme_bw() + theme(panel.grid.minor.x = element_blank(), panel.grid.major.x = element_blank())
                      + theme(panel.grid.minor.y = element_blank(), panel.grid.major.y = element_blank())
                      + theme(axis.text.x = element_text(angle = 90, hjust = 1,size=7))
                      + ggtitle(title))

    if(plotIndex > 1 || dropSigLegend)
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

plot_ratio_range_distribution<-function(ratioRangeData, bgId, freq, minRatioPerc = 0.05, minRatioSeg = 0)
{
  colCount = ncol(ratioRangeData)
  freqColStart = colCount - 100 + 1
  weightColEnd = freqColStart -1
  weightColStart = weightColEnd - 100 + 1

  bgRRData = ratioRangeData %>% filter(BgId==bgId)
  maxRatio = max(bgRRData$SigRatio)
  bgRRData = bgRRData %>% filter(SigRatio >= maxRatio * minRatioPerc)

  freqData = bgRRData[,freqColStart:colCount]
  weightData = bgRRData[,weightColStart:weightColEnd]
  xAxisValues = colnames(weightData)
  # print(xAxisValues)
  ratioSegmentNames = stri_replace_all_fixed(xAxisValues,'SW','')
  ratioSegments = as.numeric(ratioSegmentNames)

  # print(ratioSegmentNames)

  if(freq)
  {
    rrData = freqData
  }
  else
  {
    rrData = weightData
  }

  colnames(rrData) = ratioSegments

  rrDataTrans = as.data.frame(t(rrData))
  colnames(rrDataTrans) = bgRRData$Bucket
  rrDataTrans = cbind(ratioSegments, rrDataTrans)
  rownames(rrDataTrans) = NULL

  gatherIndex = ncol(rrDataTrans)
  rrDataTrans2 = gather(rrDataTrans, "Bucket", "Weight", 2:gatherIndex)
  colnames(rrDataTrans2) = c("RatioSegment", "Bucket", "Weight")

  rrDataTrans2 = rrDataTrans2 %>% filter(RatioSegment >= minRatioSeg)

  if(freq)
    label = "Count"
  else
    label = "Weight"

  maxPlotsPerPage = 20

  bucketList = rrDataTrans2 %>% group_by(Bucket) %>% count() %>% arrange(Bucket)
  bucketCount = nrow(bucketList)
  bucketPlotsPerPage = 20

  startBucketIndex = 1
  pageCount = 1

  while(startBucketIndex <= bucketCount)
  {
    endBucketIndex = startBucketIndex + bucketPlotsPerPage - 1
    bucketPlotList = bucketList[startBucketIndex:endBucketIndex, ]

    print(paste("pageCount=", pageCount, ", bucketCount=", bucketCount, ", start=", startBucketIndex, ", end=", endBucketIndex, sep=''))

    plotData = rrDataTrans2 %>% filter(Bucket %in% bucketPlotList$Bucket)
    maxWeight = max(plotData$Weight)

    rrPlot = (ggplot(data = plotData, aes(x=RatioSegment, y=Weight))
              + geom_line(aes(group=Bucket, colour=Bucket))
              + facet_wrap( ~ Bucket, ncol=3)
              + coord_cartesian(ylim = c(0, maxWeight * 1.25))
              + theme(axis.text.x = element_text(angle = 90, hjust = 1))
              + labs(title=paste("Group ", bgId, ": sample Ratios by ", label, sep=''), x="Ratio Segment", y=label, fill="Bucket"))

    print(rrPlot)

    startBucketIndex = endBucketIndex + 1
    pageCount = pageCount + 1
  }

  # return (rrPlot)
}


calc_lsqnonneg<-function(C, d)
{
  # stopifnot(is.numeric(C), is.numeric(d))
  if (!is.matrix(C) || !is.vector(d))
    stop("Argument 'C' must be a matrix, 'd' a vector.")
  m <- nrow(C); n <- ncol(C)
  if (m != length(d))
    stop("Arguments 'C' and 'd' have nonconformable dimensions.")

  tol = 10 * eps() * norm(C, type = "2") * (max(n, m) + 1)

  x  <- rep(0, n)             # initial point
  P  <- logical(n); Z <- !P   # non-active / active columns

  resid <- d - C %*% x
  w <- t(C) %*% resid
  wz <- numeric(n)

  # iteration parameters
  outeriter <- 0; it <- 0
  itmax <- 5 * n; exitflag <- 1

  while (any(Z) && any(w[Z] > tol))
  {
    outeriter <- outeriter + 1
    z <- numeric(n)
    wz <- rep(-Inf, n)
    wz[Z] <- w[Z]
    im <- which.max(wz)
    P[im] <- TRUE; Z[im] <- FALSE
    z[P] <- qr.solve(C[, P], d)

    while (any(z[P] <= 0))
    {
      it <- it + 1
      if(it > itmax)
      {
        print(paste("iteration count exceed at i=", i, sep=''))
        break;
        # stop("Iteration count exceeded")
      }

      Q <- (z <= 0) & P
      alpha <- min(x[Q] / (x[Q] - z[Q]))
      x <- x + alpha*(z - x)
      Z <- ((abs(x) < tol) & P) | Z
      P <- !Z
      z <- numeric(n)
      z[P] <- qr.solve(C[, P], d)
    }
    x <- z
    resid <- d - C %*% x
    w <- t(C) %*% resid
  }
  return(list(x = x, resid.norm = sum(resid*resid)))
}
