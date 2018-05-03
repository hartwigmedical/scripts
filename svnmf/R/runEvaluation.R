evaluate_run<-function(runNumber, sigCount, nmfMatrixData, summaryCounts, signatures, contribution, bucketNames, sigNames, runEstimates = FALSE) {

  # prepare PDF for output
  pdf(file=paste("~/logs/r_output/svnmf_", runNumber, ".pdf", sep = ""), height = 14, width = 20)

  sampleIds = summaryCounts[,1]

  if(runEstimates == TRUE) {
    # estimate functions
    nmfEstimate <- nmf(nmfMatrixData, rank=6:12, method="brunet", nrun=4, seed=123456)
    plot(nmfEstimate)
  }

  # Bucket Evaluation
  sigBucketData = svnmf::get_bucket_data(signatures, contribution, bucketNames)

  sigBucketStats = svnmf::get_sig_bucket_stats(sigBucketData)

  # Signature Discovery, by looking at relative contribution of buckets

  # report on top 3 contributing buckets only
  sigBucketTop3 = svnmf::get_top_buckets(sigBucketData, 3)

  # key bucket stats
  bucketSummaryData = svnmf::get_bucket_stats(sigBucketData)

  # least contributing 10 buckets
  leastContribBuckets = svnmf::get_least_contrib_buckets(bucketSummaryData)

  # output all bucket-related analysis
  ggarrange(ggtexttable(sigBucketStats, rows = NULL),
            ggtexttable(sigBucketTop3, rows = NULL),
            ggtexttable(bucketSummaryData, rows = NULL),
            ggtexttable(leastContribBuckets, rows = NULL),
            ncol = 1, nrow = 4, heights = c(1, 1, 1, 1))


  # compare to signature graphs

  # plot as a bar chart - first put the data back into the required format and add bucket names
  df = data.frame(bucket = bucketNames)
  df2 = cbind(df, sigPercents)
  sigPlotData = melt(df2, id.vars = c("bucket"))

  ymax = 1.0

  sigBucketsPlot = (ggplot(data = sigPlotData, aes(x = bucket, y = value, fill = bucket))
                  + geom_bar(stat = "identity", colour = "black", size = 0.2)
                  + facet_wrap(~ variable,ncol = 1) + ylab("Relative contribution")
                  + coord_cartesian(ylim = c(0, ymax))
                  + theme(axis.text.x = element_text(angle = 90, hjust = 1))
                  + scale_y_continuous(breaks = seq(0, ymax, 0.1))
                  + ggtitle("Signatures by Bucket")
                  + guides(fill = FALSE))

  ggarrange(sigBucketsPlot, ncol = 1, nrow = 1)

  # 2 Signature Evaluation
  sampleSigData = svnmf::get_sig_data(signatures, contribution, sigNames, sampleNames)

  # key stats per signature
  sigStats = svnmf::get_sig_stats(sampleSigData)


  ggarrange(ggtexttable(sigStats, rows = NULL), ncol = 1, nrow = 1)

  # run again, this time bucketing samples into mutational load



  # grid.table(sigStats, rows = NULL)

  dev.off()



  # other NMF output


}
