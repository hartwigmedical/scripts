evaluate_nmf_run<-function(runType, runId, sigCount, nmfResult, matrixData, summaryCounts, sampleCancerTypes, bucketNames,
                           sigNamesUnamed, sigNamesNamed, plotByCancerType = T, viewResults = F, calcResiduals = F)
{
  # use this version when the NMF result is not available (eg when using external signatures like COSMIC)
  signatures = NMF::basis(nmfResult)
  contribution = NMF::coef(nmfResult)

  evaluate_nmf_data(runType, runId, sigCount, signatures, contribution, matrixData, summaryCounts,
                    sampleCancerTypes, bucketNames, sigNamesUnamed, sigNamesNamed, plotByCancerType, viewResults, calcResiduals)
}

evaluate_nmf_data<-function(runType, runId, sigCount, signatures, contribution, matrixData, summaryCounts, sampleCancerTypes, bucketNames,
                            sigNamesUnamed, sigNamesNamed, plotByCancerType = T, viewResults = F, calcResiduals = T)
{
  print(paste("evaluating run: type=", runType, ", id=", runId, ", sigCount=", sigCount, sep=''))

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
  # leastContribBuckets = get_least_contrib_buckets(bucketSummaryData)
  # leastContribBuckets = head(leastContribBuckets, 10)

  if(viewResults)
  {
    # View(sigBucketData)
    View(sigBucketStats)
    View(sigBucketTopN)
  }

  # Top Bucket Counts per Sample - for now get all buckets
  # sampleBucketTopN = get_top_buckets_by_sample(summaryCounts, origSampleCounts, sampleCancerTypes, 0)
  sampleBucketData = get_sample_bucket_data(matrixData, origSampleCounts, bucketNames)

  # 2 Signature Evaluation
  print("evaluating signatures")

  # optionally name signatues for subsequent output
  sigNames = sigNamesNamed

  sampleSigData = get_sig_data(signatures, contribution, sigNames, sampleNames)

  # key stats per signature
  sigStats = get_sig_stats(sampleSigData)

  # calculate and factor in residuals
  if(calcResiduals)
  {
    sampleSigData = append_residuals(contribution, signatures, matrixData, bucketNames, sampleSigData)
  }

  # run again, this time bucketing samples into mutational load and cancer types

  # get cancer type and SV Count
  sampleSigData = merge(sampleSigData, sampleCancerTypes,by.x="SampleId",by.y="SampleId",all.x=TRUE)
  sampleSigData$CancerType = ifelse(is.na(sampleSigData$CancerType), 'N/A', paste(sampleSigData$CancerType, sep=""))

  sampleSigCounts = sampleSigData %>% filter(SigName!="Residual") %>% group_by(SampleId) %>% summarise(SampleCount=sum(Count))
  sampleSigData = merge(sampleSigData, sampleSigCounts, by.x="SampleId",by.y="SampleId",all.x=TRUE)

  sampleSigDataNoRes = sampleSigData %>% filter(SigName!="Residual")

  if(viewResults)
  {
    View(sigStats)
    View(sampleSigData)
  }

  outputFile = paste("~/logs/r_output/pdfs/", runType, "_", runId, ".pdf", sep = "")
  print(paste("writing output to file: ", outputFile, sep=''))

  # DATA OUTPUT TO PDF
  pdf(file=outputFile, height = 14, width = 20)

  par(mar=c(1,1,1,1))

  # 1. Bucket data
  title = textGrob("Bucket Summary Data & Top-N Buckets", gp=gpar(fontface="bold", fontsize=16))
  plot_bucket_summary_data(bucketSummaryData, sigBucketTopN, title)

  sigColours = get_sig_colours(sigCount)

  sigColoursPlusRes = sigColours
  sigColoursPlusRes[sigCount+1] = "#000000"

  bucketSummaryPlot = get_bucket_summary_plot(bucketSummaryData, runType)
  grid.arrange(bucketSummaryPlot, ncol = 1, nrow = 1, newpage = TRUE)

  title = textGrob("Signature-Bucket Stats", gp=gpar(fontface="bold", fontsize=16)) #  & Least Important Buckets
  grid.arrange(tableGrob(sigBucketStats, rows=NULL),
               # tableGrob(leastContribBuckets, rows = NULL),
               widths = c(2, 1), ncol = 2, newpage = TRUE, top=title)

  # 2. Default signature-bucket plot
  sigBucketsPlots = get_bucket_signatures_plot(bucketNames, signatures, sigNames)
  for(i in 1:length(sigBucketsPlots))
  {
    grid.arrange(sigBucketsPlots[[i]], ncol=1, nrow=1, newpage = TRUE)
  }

  # 3. Plot of Bucket Contribution for each Sample
  plot_sample_bucket_contrib(sampleBucketData, nrow(bucketNames), runType, 0) # sampleBucketTopN

  # 4. Signature data
  title = textGrob("Signature Stats", gp=gpar(fontface="bold", fontsize=16))
  grid.arrange(tableGrob(sigStats, rows = NULL), ncol = 1, nrow = 1, top=title, newpage = TRUE)

  # 5. Relative Signature contributions by Cancer Type
  plot_cancer_sigs(sampleSigDataNoRes, sigColours, runType)

  # 6. Top 50 samples by signature, but include all other signatures as well
  plot_top_n_samples_by_sig(sampleSigDataNoRes, sigNames, 50, sigColours, runType)

  # 7. Sigs with Samples by cancer type
  plot_sig_samples(sampleSigData, "", sigColoursPlusRes, runType) # all samples

  cancerTypes = sampleCancerTypes %>% group_by(CancerType) %>% count()

  if(plotByCancerType)
  {
    for(cancerType in cancerTypes$CancerType)
    {
      if(!is.na(cancerType))
      {
        plot_sig_samples(sampleSigData, cancerType, sigColoursPlusRes, runType)
      }
    }
  }

  dev.off()
}




