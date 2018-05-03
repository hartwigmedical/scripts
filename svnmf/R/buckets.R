get_bucket_data<-function(signatures, contribution, bucketNames) {

  # work out percentage for each bucket's contribution
  sigPercents = apply(signatures, 2, function(x) x/sum(x))
  sigPercentsByBuckets = melt(cbind(data.frame(bucket = bucketNames), sigPercents), id.vars = c("bucket"))
  colnames(sigPercentsByBuckets) <- c('Bucket', 'Signature', 'Percent')

  # use signature counts to turn contributions into counts
  sigBucketCounts = signatures
  for (i in 1:ncol(sigBucketCounts))
  {
    sigBucketCounts[,i] = round(sigBucketCounts[,i]* rowSums(contribution)[i],0)
  }

  sigByBuckets = melt(cbind(data.frame(bucket = bucketNames), sigBucketCounts), id.vars = c("bucket"))
  colnames(sigByBuckets) <- c('Bucket', 'Signature', 'SvCount')

  # merge counts and percentages into a single signature by bucket dataframe
  sigBucketData = cbind(sigPercentsByBuckets, SvCount = round(sigByBuckets$SvCount,0))
  return(sigBucketData)
}

get_sig_bucket_stats<-function(sigBucketData) {

  sigBucketStats = (sigBucketData %>% group_by(Signature)
                    %>% summarise(BucketCount=sum(SvCount>0),
                                  BucketPerc=round(sum(SvCount>0)/n_distinct(Bucket),2),
                                  SvCount=sum(SvCount),
                                  MaxPercent=round(max(Percent),3),
                                  AvgPercent=round(sum(ifelse(Percent>0.001,Percent,0))/sum(Percent>0),3),
                                  PercGT75=sum(Percent>0.75),
                                  Perc50_75=sum(Percent>0.5&Percent<=0.75),
                                  Perc25_50=sum(Percent>0.25&Percent<=0.5),
                                  Perc10_25=sum(Percent>0.1&Percent<=0.25),
                                  Perc5_10=sum(Percent>0.05&Percent<=0.1),
                                  PercLT5=sum(Percent>0.001&Percent<=0.05))
                    %>% arrange(-BucketCount))

  return (sigBucketStats)
}

get_top_buckets<-function(sigBucketData, sigNames, sigNamesNamed, maxN = 4, reqPercent = 0.85) {

  sigBucketTopN = data.frame()
  for(i in 1:length(sigNames))
  {
    # take top 3 entries for each signature in turn, and append them to the result
    sigData = sigBucketData %>% filter(Signature==i) %>% arrange(-Percent) # sig-agnostic version

    percentTotal = 0
    for(i in 1:nrow(sigData))
    {
      bucketRow = sigData[i,]
      sigBucketTopN = rbind(sigBucketTopN, bucketRow)

      percentTotal = percentTotal + bucketRow$Percent
      if(percentTotal >= reqPercent | i >= maxN)
      {
        break
      }
    }
  }

  sigBucketTopN$Percent = round(sigBucketTopN$Percent, 3)

  # meld on named sigs
  sigNamesCombined = cbind(sigNames, sigNamesNamed)
  colnames(sigNamesCombined) <- c("Signature", "SigName")
  sigBucketTopN = (merge(sigBucketTopN, sigNamesCombined, by.x="Signature", by.y="Signature", all.x=TRUE) %>% arrange(SigName,-Percent))
  return (sigBucketTopN)
}

get_bucket_stats<-function(sigBucketData) {
  bucketSummaryData = (sigBucketData %>% group_by(Bucket)
                       %>% summarise(SigCount=sum(SvCount>0),
                                     SigPerc=round(sum(SvCount>0)/n_distinct(Signature),2),
                                     SvCount=sum(SvCount),
                                     MaxPercent=round(max(Percent),3),
                                     AvgPercent=round(sum(ifelse(Percent>0.001,Percent,0))/sum(Percent>0),3))
                       %>% arrange(-SvCount))

  return (bucketSummaryData)
}

get_least_contrib_buckets<-function(bucketSummaryData, worstN = 20) {
  return (top_n(bucketSummaryData, worstN, -SvCount) %>% arrange(SvCount))
}

get_bucket_signatures_plot<-function(bucketNames, signatures, sigNames) {

  df = data.frame(Bucket = bucketNames)
  sigPercents = apply(signatures, 2, function(x) x/sum(x))
  colnames(sigPercents) <- sigNames
  df2 = cbind(df, sigPercents)
  sigPlotData = melt(df2, id.vars = c("Bucket"))

  ymax = 1.0

  sigBucketsPlot = (ggplot(data = sigPlotData, aes(x = Bucket, y = value, fill = Bucket))
                    + geom_bar(stat = "identity", colour = "black", size = 0.2)
                    + facet_wrap(~ variable,ncol = 1)
                    + ylab("Relative contribution")
                    + coord_cartesian(ylim = c(0, ymax))
                    + theme(axis.text.x = element_text(angle = 90, hjust = 1))
                    + ggtitle("Signatures by Bucket")
                    + guides(fill = FALSE))

  return (sigBucketsPlot)
}

get_bucket_summary_plot<-function(bucketSummaryData) {

    bucketSummaryPlot <- (ggplot(data = bucketSummaryData,
                               aes(x=reorder(Bucket, -SvCount), y = SvCount, group = 1), fill = Bucket)
                        + geom_bar(stat = "identity", colour = "black", size = 0.2)
                        + ylab("SV Count")
                        + xlab("Bucket")
                        + theme(axis.text.x = element_text(angle = 90, hjust = 1))
                        + ggtitle("Bucket Summary"))

  axisRatio = max(bucketSummaryData$SvCount) / max(bucketSummaryData$SigCount) * 0.6

  bucketSummaryPlot <- (bucketSummaryPlot + geom_line(aes(y = SigCount*axisRatio, color = "red"))
                        + scale_y_continuous(sec.axis = sec_axis(~.*(1/axisRatio), name = "Sig Count"))
                        + theme(legend.position="none"))

  return (bucketSummaryPlot)
}

