get_bucket_data<-function(signatures, contribution, bucketNames)
{
  # work out percentage for each bucket's contribution
  sigPercents = apply(signatures, 2, function(x) x/sum(x))

  sigPercents = cbind(bucketNames,sigPercents)
  sigPercentsByBuckets = gather(as.data.frame(sigPercents), "Signature", "Percent", -Bucket)

  # use signature counts to turn contributions into counts
  sigBucketCounts = signatures
  for (i in 1:ncol(sigBucketCounts))
  {
    sigBucketCounts[,i] = round(sigBucketCounts[,i]* rowSums(contribution)[i],0)
  }

  sigBucketCounts = cbind(bucketNames,sigBucketCounts)
  sigCountsByBuckets = gather(as.data.frame(sigBucketCounts), "Signature", "Count", -Bucket)

  # finally merge the sig percents and counts by bucket
  sigBucketData = merge(sigPercentsByBuckets,sigCountsByBuckets,by=c("Signature","Bucket"),all.x=T)
  return(sigBucketData)
}

get_sig_bucket_stats<-function(sigBucketData) {

  sigBucketStats = (sigBucketData %>% group_by(SigName)
                    %>% summarise(BucketCount=sum(Count>0),
                                  BucketPerc=round(sum(Count>0)/n_distinct(Bucket),2),
                                  Count=sum(Count),
                                  MaxPercent=round(max(Percent),3),
                                  AvgPercent=round(sum(ifelse(Percent>0.001,Percent,0))/sum(Percent>0),3),
                                  PercGT75=sum(Percent>0.75),
                                  Perc50_75=sum(Percent>0.5&Percent<=0.75),
                                  Perc25_50=sum(Percent>0.25&Percent<=0.5),
                                  Perc10_25=sum(Percent>0.1&Percent<=0.25),
                                  Perc5_10=sum(Percent>0.05&Percent<=0.1),
                                  PercLT5=sum(Percent>0.001&Percent<=0.05))
                    %>% arrange(SigName))

  return (sigBucketStats)
}

get_top_buckets<-function(sigBucketData, maxN = 4, reqPercent = 0.85)
{
  sigBucketTopN = data.frame()
  sigCount = n_distinct(sigBucketData$Signature)
  for(i in 1:sigCount)
  {
    # take top 3 entries for each signature in turn, and append them to the result
    sigData = sigBucketData %>% filter(Signature==i) %>% arrange(-Percent) # sig-agnostic version

    percentTotal = 0
    for(j in 1:nrow(sigData))
    {
      bucketRow = sigData[j,]
      sigBucketTopN = rbind(sigBucketTopN, bucketRow)

      percentTotal = percentTotal + bucketRow$Percent
      if(percentTotal >= reqPercent | j >= maxN)
      {
        break
      }
    }
  }

  sigBucketTopN$Percent = round(sigBucketTopN$Percent, 3)

  # merge on named sigs
  sigBucketTopN = sigBucketTopN %>% arrange(SigName,-Percent)
  # sigBucketTopN = merge(sigBucketTopN, sigNamesCombined, by.x="Signature", by.y="Signature", all.x=T) %>% arrange(SigName,-Percent)
  return (sigBucketTopN)
}

get_bucket_stats<-function(sigBucketData) {
  bucketSummaryData = (sigBucketData %>% group_by(Bucket)
                       %>% summarise(SigCount=sum(Count>0),
                                     SigPerc=round(sum(Count>0)/n_distinct(Signature),2),
                                     Count=sum(Count),
                                     MaxPercent=round(max(Percent),3),
                                     AvgPercent=round(sum(ifelse(Percent>0.001,Percent,0))/sum(Percent>0),3))
                       %>% arrange(-Count))

  return (bucketSummaryData)
}

get_least_contrib_buckets<-function(bucketSummaryData, worstN = 20) {
  return (top_n(bucketSummaryData, worstN, -Count) %>% arrange(Count))
}

plot_bucket_summary_data<-function(bucketSummaryData, sigBucketTopN, title, rowsPerColumn)
{
  topNRows = nrow(sigBucketTopN)
  if(topNRows <= rowsPerColumn)
  {
    grid.arrange(tableGrob(head(bucketSummaryData, rowsPerColumn), rows=NULL),
                 tableGrob(sigBucketTopN, rows=NULL),
                 ncol = 2, newpage = TRUE, top=title)
  }
  else
  {
    # print the first 40 (or whatever RPC is) on the page with the bucket stats
    grid.arrange(tableGrob(head(bucketSummaryData, rowsPerColumn), rows=NULL),
                 tableGrob(head(sigBucketTopN,rowsPerColumn), rows=NULL),
                 ncol = 2, newpage = TRUE, top=title)

    i = rowsPerColumn+1
    while(i < topNRows)
    {
      # if(i + (rowsPerColumn*2) - 1 <= topNRows)
      if(topNRows > i + rowsPerColumn-1) # say i = 41 then 2 columns required iftopNRows > i + 40-1
        {
        # if 2 columns on the new page are required
        s1 = i
        s2 = s1 + rowsPerColumn-1
        s3 = s2+1
        s4 = min(s3 + rowsPerColumn-1, topNRows)

        # print(paste(i, ", topNRows=", topNRows, ", s1=", s1, ", s2=", s2, ", s3=", s3, ", s4=", s4, sep=''))

        grid.arrange(tableGrob(sigBucketTopN[s1:s2,], rows=NULL),
                     tableGrob(sigBucketTopN[s3:s4,], rows=NULL),
                     ncol = 2, newpage = TRUE, top=title)
        i = s4+1
      }
      else
      {
        # otherwise only 1 page and column left to print
        maxRow = min(i+rowsPerColumn-1, topNRows)
        # print(paste(i, ", topNRows=", topNRows, ", maxRow=", maxRow, sep=''))
        grid.arrange(tableGrob(sigBucketTopN[i:maxRow,], rows=NULL),
                     ncol = 2, newpage = TRUE, top=title)

        i = topNRows
        break
      }
    }
  }

}

get_bucket_signatures_plot<-function(bucketNames, signatures, sigNames) {

  df = data.frame(Bucket = bucketNames)
  sigPercents = apply(signatures, 2, function(x) x/sum(x))
  colnames(sigPercents) <- sigNames
  df2 = cbind(df, sigPercents)
  sigPlotData = melt(df2, id.vars = c("Bucket"))

  ymax = min(1.0, max(sigPlotData$value)+0.1)
  bucketCount = nrow(bucketNames)
  dataCount = nrow(sigPlotData)
  sigCount = ncol(signatures)
  bucketsPerPage = 10
  rowsPerPage = bucketCount * bucketsPerPage
  rowStart = 1

  sigBucketPlots = list()
  plotIndex = 1

  while(rowStart < dataCount)
  {
    if(sigCount >= 15)
    {
      rowEnd = min(rowStart + rowsPerPage - 1, dataCount)
    }
    else
    {
      rowEnd = dataCount
    }

    sigBucketPlot = (ggplot(data = sigPlotData[rowStart:rowEnd,], aes(x = Bucket, y = value, fill = Bucket))
                      + geom_bar(stat = "identity", colour = "black", size = 0.2)
                      + facet_wrap(~ variable,ncol = 1)
                      + ylab("Relative contribution")
                      + coord_cartesian(ylim = c(0, ymax))
                      + theme(axis.text.x = element_text(angle = 90, hjust = 1))
                      + ggtitle("Signatures by Bucket")
                      + guides(fill = FALSE))

    sigBucketPlots[[plotIndex]] = sigBucketPlot
    plotIndex = plotIndex + 1

    rowStart = rowEnd + 1
  }


  return (sigBucketPlots)
}

get_bucket_summary_plot<-function(bucketSummaryData, varType = "SV") {

    bucketSummaryPlot <- (ggplot(data = bucketSummaryData,
                               aes(x=reorder(Bucket, -Count), y = Count, group = 1), fill = Bucket)
                        + geom_bar(stat = "identity", colour = "black", size = 0.2)
                        + ylab(paste(varType, " Count", sep=''))
                        + xlab("Bucket")
                        + theme(axis.text.x = element_text(angle = 90, hjust = 1))
                        + ggtitle("Bucket Summary"))

  axisRatio = max(bucketSummaryData$Count) / max(bucketSummaryData$SigCount) * 0.6

  bucketSummaryPlot <- (bucketSummaryPlot + geom_line(aes(y = SigCount*axisRatio, color = "red"))
                        + scale_y_continuous(sec.axis = sec_axis(~.*(1/axisRatio), name = "Sig Count"))
                        + theme(legend.position="none"))

  return (bucketSummaryPlot)
}

get_sample_bucket_data<-function(matrixData, origSampleCounts, bucketNames)
{
  # convert the matrix back into a dataframe organise by SampleId and Bucket
  # and ensure that every bucket is in every sample, even if the count is zero
  data1 = t(matrixData)
  data2 = t(data1)
  data3 = cbind(Bucket = bucketNames, data2)
  data4 = gather(data3, SampleId, Count, -Bucket)

  data4 = merge(data4, origSampleCounts, by.x="SampleId",by.y="SampleId",all.x=TRUE)
  colnames(data4) <- c("SampleId", "Bucket", "Count", "SampleCount")
  data4$Percent = round(data4$Count/data4$SampleCount,4)

  data4[is.na(data4)] = 0

  return (data4)
}

get_top_buckets_by_sample<-function(sampleCounts, origSampleCounts, sampleCancerTypes, topNBuckets = 20)
{
  sampleBucketData = merge(sampleCounts, origSampleCounts, by.x="SampleId",by.y="SampleId",all.x=TRUE)
  sampleBucketData = merge(sampleBucketData, sampleCancerTypes,by.x="SampleId",by.y="SampleId",all.x=TRUE)
  sampleBucketData = setNames(sampleBucketData, c("SampleId", "Bucket", "Count", "SampleCount", "CancerType"))

  sampleBucketData$Percent = round(sampleBucketData$Count/sampleBucketData$SampleCount,4)

  # collect the top N buckets by count for each sample
  sampleBucketData = sampleBucketData %>% arrange(SampleId,-Count)

  if(topNBuckets == 0)
  {
    return (sampleBucketData)
  }

  sampleBucketTopN = data.frame()
  for(sampleId in sampleCancerTypes$SampleId)
  {
    samBucketData = head(sampleBucketData %>% filter(SampleId==sampleId),topNBuckets)
    sampleBucketTopN = rbind(sampleBucketTopN,samBucketData)
  }

  return (sampleBucketTopN)
}

plot_sample_bucket_contrib<-function(sampleBucketData, bucketCount, varType = "SV", maxPlots = 4, extTitle = "", isAbsolute = F)
{
  bucketColours = c("#ff994b", "#463ec0", "#88c928", "#996ffb", "#68b1c0", "#e34bd9", "#106b00", "#d10073", "#98d76a", "#6b3a9d",
                   "#d5c94e", "#0072e2", "#ff862c", "#31528d", "#d7003a", "#323233", "#ff4791", "#01837a", "#ff748a", "#777700",
                   "#ff86be", "#4a5822", "#ffabe4", "#6a4e03", "#c6c0fb", "#ffb571", "#873659", "#dea185", "#a0729d", "#8a392f",
                   "#ff984b", "#469ec0", "#88c926", "#997ffb", "#67b1c0", "#e35bd9", "#105b00", "#d11073", "#98676a", "#6b9a9d",
                   "#d5c84e", "#0092e2", "#ff8626", "#31728d", "#d6003a", "#325233", "#ff5791", "#01137a", "#ff648a", "#779700",
                   "#ff88be", "#4a9822", "#ffabe6", "#6a7e03", "#c6c0fb", "#ff5571", "#875659", "#de1185", "#a0529d", "#8a992f",
                   "#ff974b", "#468ec0", "#88c925", "#998ffb", "#65b1c0", "#e36bd9", "#107b00", "#d12073", "#98576a", "#6b8a9d",
                   "#d5c74e", "#0082e2", "#ff8625", "#31828d", "#d5003a", "#326233", "#ff7791", "#01237a", "#ff548a", "#778700",
                   "#ff87be", "#4a8822", "#ffabe5", "#6a8e03", "#c5c0fb", "#ff6571", "#877659", "#de2185", "#a0529d", "#8a892f",
                   "#ff964b", "#467ec0", "#88c924", "#999ffb", "#64b1c0", "#e37bd9", "#108b00", "#d13073", "#98476a", "#6b4a9d",
                   "#d5c64e", "#0062e2", "#ff8624", "#31928d", "#d4003a", "#327233", "#ff8791", "#01337a", "#ff448a", "#774700",
                   "#ff86be", "#4a7822", "#ffabe4", "#6a9e03", "#c4c0fb", "#ff7571", "#878659", "#de3185", "#a0429d", "#8a492f")


  sbData = sampleBucketData %>% filter(SampleCount>0)

  sbData = sbData %>% arrange(-SampleCount, SampleId) %>% select('SampleId', 'Bucket', 'Percent', "Count", "SampleCount")

  bucketPlots = list()
  plotIndex = 1

  if(nrow(sbData) > 0)
  {
    # only plot 50 samples at a time
    samplesPerPlot = 100
    rowsPerPlot = samplesPerPlot * bucketCount

    rowEnd = 0
    maxRows = nrow(sbData)
    plotCount = 0

    while(rowEnd < maxRows)
    {
      if(plotCount == 0)
      {
        rowStart = 1
        rowEnd = min(rowStart + 80*bucketCount - 1, maxRows) # less to leave space for the legend
      }
      else
      {
        rowStart = rowEnd + 1
        rowEnd = min(rowStart + rowsPerPlot - 1, maxRows)

        # combine stragglers onto the last graph
        if(rowEnd < maxRows & maxRows - rowEnd < rowsPerPlot * 0.4)
        {
          # print(paste(plotCount, ": rowEnd=", rowEnd, " close to maxRows=", maxRows, sep=''))
          rowEnd = maxRows
        }
      }

      if(extTitle == "")
        title = "Bucket % by Sample"
      else
        title = extTitle

      if(isAbsolute)
      {
        sampleBucketPlot <- ggplot(sbData[rowStart:rowEnd,], aes(x = reorder(SampleId, -SampleCount), y = Count, fill = Bucket))
      }
      else
      {
        sampleBucketPlot <- ggplot(sbData[rowStart:rowEnd,], aes(x = reorder(SampleId, -SampleCount), y = Percent, fill = Bucket))
      }

      sampleBucketPlot <- (sampleBucketPlot
                         + geom_bar(stat = "identity", colour = "black", size=0.01)
                         + geom_bar(stat = "identity", colour = "black", size=0.01)
                         + labs(x = "", y = paste("Bucket % by Sample", sep=''))
                         + scale_fill_manual(values = bucketColours)
                         + theme_bw() + theme(panel.grid.minor.x = element_blank(), panel.grid.major.x = element_blank())
                         + theme(panel.grid.minor.y = element_blank(), panel.grid.major.y = element_blank())
                         + theme(axis.text.x = element_text(angle = 90, hjust = 1,size=7)))

      if(plotCount == 0)
      {
        sampleBucketPlot <- sampleBucketPlot + ggtitle(title)
      }
      else
      {
        # remove legend after the first plot
        sampleBucketPlot <- sampleBucketPlot + theme(legend.position="none")
      }

      bucketPlots[[plotIndex]] <- sampleBucketPlot
      plotCount = plotCount + 1

      if(plotIndex >= 2)
      {
        multiplot(plotlist = bucketPlots, cols = 1)
        bucketPlots = list()
        plotIndex = 1
      }
      else
      {
        plotIndex = plotIndex + 1
      }

      if(maxPlots > 0 && plotCount >= maxPlots)
        break
    }

    if(plotIndex > 1)
    {
      # now print all plots for this cancer type
      multiplot(plotlist = bucketPlots, cols = 1)
    }
  }
}

