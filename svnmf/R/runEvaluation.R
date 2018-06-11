evaluate_nmf_run<-function(runType, runId, sigCount, nmfResult, summaryCounts, sampleCancerTypes, bucketNames,
                           sigNamesUnamed, sigNamesNamed, plotByCancerType = T, viewResults = F)
{
  print(paste("evaluating run: type=", runType, ", id=", runId, ", sigCount=", sigCount, sep=''))

  signatures = NMF::basis(nmfResult)
  contribution = NMF::coef(nmfResult)
  sampleNames = colnames(contribution)

  origSampleCounts = summaryCounts %>% group_by(SampleId) %>% summarise(OrigSampleCount=sum(Count))

  print("evaluating buckets")

  # 1. Bucket Evaluation
  sigBucketData = get_bucket_data(signatures, contribution, bucketNames)
  sigBucketStats = get_sig_bucket_stats(sigBucketData)

  # Signature Discovery, by looking at relative contribution of buckets
  sigBucketTopN = get_top_buckets(sigBucketData, sigNamesUnamed, sigNamesNamed)

  # key bucket stats
  bucketSummaryData = get_bucket_stats(sigBucketData)

  # least contributing 10 buckets
  leastContribBuckets = get_least_contrib_buckets(bucketSummaryData)
  leastContribBuckets = head(leastContribBuckets, 10)

  if(viewResults)
  {
    # View(sigBucketData)
    View(sigBucketStats)
    View(sigBucketTopN)
  }

  # Top Bucket Counts per Sample
  sampleBucketTopN = get_top_buckets_by_sample(summaryCounts, origSampleCounts, sampleCancerTypes)

  # 2 Signature Evaluation
  print("evaluating signatures")

  # optionally name signatues for subsequent output
  sigNames = sigNamesNamed

  sampleSigData = get_sig_data(signatures, contribution, sigNames, sampleNames)

  # key stats per signature
  sigStats = get_sig_stats(sampleSigData)

  # run again, this time bucketing samples into mutational load and cancer types

  # get cancer type and SV Count
  sampleSigData = merge(sampleSigData, sampleCancerTypes,by.x="SampleId",by.y="SampleId",all.x=TRUE)
  sampleSigData$CancerType = ifelse(is.na(sampleSigData$CancerType), 'N/A', paste(sampleSigData$CancerType, sep=""))

  sampleSigCounts = sampleSigData %>% group_by(SampleId) %>% summarise(SampleCount=sum(Count))
  sampleSigData = merge(sampleSigData, sampleSigCounts, by.x="SampleId",by.y="SampleId",all.x=TRUE)

  if(viewResults)
  {
    View(sigStats)
    View(sampleSigData)
  }

  outputFile = paste("~/logs/r_output/", runType, "_", runId, ".pdf", sep = "")
  print(paste("writing output to file: ", outputFile, sep=''))

  # DATA OUTPUT TO PDF
  pdf(file=outputFile, height = 14, width = 20)

  par(mar=c(1,1,1,1))

  # 2. bucket data
  title = textGrob("Bucket Summary Data & Top-N Buckets", gp=gpar(fontface="bold", fontsize=16))
  plot_bucket_summary_data(bucketSummaryData, sigBucketTopN)

  sigColours = get_sig_colours(sigCount)

  bucketSummaryPlot = get_bucket_summary_plot(bucketSummaryData, runType)
  grid.arrange(bucketSummaryPlot, ncol = 1, nrow = 1, newpage = TRUE)

  title = textGrob("Signature-Bucket Stats & Least Important Buckets", gp=gpar(fontface="bold", fontsize=16))
  grid.arrange(tableGrob(sigBucketStats, rows=NULL),
               tableGrob(leastContribBuckets, rows = NULL),
               ncol = 2, newpage = TRUE, top=title)

  # 3. default signature-bucket plot
  sigBucketsPlot = get_bucket_signatures_plot(bucketNames, signatures, sigNames)
  grid.arrange(sigBucketsPlot, ncol = 1, nrow = 1, newpage = TRUE)

  # Plot of Bucket Contribution for each Sample
  plot_sample_bucket_contrib(sampleBucketTopN, "", 10, runType, 4)

  # 4. Signature data
  title = textGrob("Signature Stats", gp=gpar(fontface="bold", fontsize=16))
  grid.arrange(tableGrob(sigStats, rows = NULL), ncol = 1, nrow = 1, top=title, newpage = TRUE)

  # 5. Relative Signature contributions by Cancer Type
  plot_cancer_sigs(sampleSigData, sigColours, runType)

  # 6. Top 50 samples by signature, but include all other signatures as well
  plot_top_n_samples_by_sig(sampleSigData, sigNames, 50, sigColours, runType)

  # 7. Sigs with Samples by cancer type
  plot_sig_samples(sampleSigData, "", sigColours, runType) # all samples

  cancerTypes = sampleCancerTypes %>% group_by(CancerType) %>% count()

  if(plotByCancerType)
  {
    for(cancerType in cancerTypes$CancerType)
    {
      if(!is.na(cancerType))
      {
        plot_sig_samples(sampleSigData, cancerType, sigColours, runType)
      }
    }
  }

  dev.off()
}
