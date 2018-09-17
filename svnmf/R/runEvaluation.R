evaluate_nmf_run<-function(runType, runId, nmfResult, matrixData, summaryCounts, sampleCancerTypes, bucketNames,
                           sigNamesNamed, plotByCancerType = T, viewResults = F)
{
  # use this version when the NMF result is not available (eg when using external signatures like COSMIC)
  signatures = NMF::basis(nmfResult)
  contribution = NMF::coef(nmfResult)

  evaluate_nmf_data(runType, runId, signatures, contribution, matrixData, summaryCounts,
                    sampleCancerTypes, bucketNames, sigNamesNamed, plotByCancerType, viewResults)
}

evaluate_nmf_data<-function(runType, runId, signatures, contribution, matrixData, summaryCounts, sampleCancerTypes, bucketNames,
                            sigNamesNamed, plotByCancerType = T, viewResults = F, bgSigCount = 0, printAllPlots = T)
{
  evaluate_signature_fit(runType, runId, signatures, contribution, matrixData, summaryCounts, sampleCancerTypes, bucketNames, sigNamesNamed, plotByCancerType, viewResults)
}

produce_signature_report<-function(runType, runId, sigFile, contribFile, sigInfoFile,
                                   matrixData, summaryCounts, sampleCancerTypes, bucketNames,
                                   plotByCancerType = T, viewResults = F, printAllPlots = T)
{
  baContribs = as.matrix(read.csv(file=contribFile, stringsAsFactors=F))
  baSigs = as.matrix(read.csv(file=sigFile, stringsAsFactors=F))
  baSigNames = colnames(baSigs)
  baSigCount = ncol(baSigs)

  baSigNumList = get_signame_list(baSigCount, F)
  baSigStrList = get_signame_list(baSigCount, T)
  colnames(baSigs) <- baSigNumList

  # trim sig names
  for(i in 1:baSigCount)
  {
    baSigNames[i] = paste(baSigStrList[i], baSigNames[i], sep="_")
    sigLen = stringi::stri_length(baSigNames[i])

    if(sigLen > 8)
    {
      catIndex = stri_locate_first_fixed(baSigNames[i], "_cat")
      if(!is.na(catIndex[1]))
      {
        baSigNames[i] = substring(baSigNames[i], 1, catIndex[1]-1)
      }

      baSigNames[i] = stri_replace_all_fixed(baSigNames[i], ".", "")
    }
  }

  baSigInfo = as.data.frame(read.csv(file=sigInfoFile, stringsAsFactors=F))

  # look for links between the sigs
  baSigCount = nrow(baSigInfo)

  baseColours = get_base_colours()

  # look for links between the sigs
  bgSigCount = 0
  majorSigCount = 0
  newSigColours = sigColours

  for(i in 1:baSigCount)
  {
    sigInfo = baSigInfo[i,]

    if(sigInfo$Type == "BGRD")
    {
      bgSigCount = bgSigCount + 1
      newSigColours[i] = "grey60"
    }
    else if(sigInfo$Type == "MAJOR" && (is.na(sigInfo$ParentId) || sigInfo$ParentId==""))
    {
      bgId = sigInfo$BgId
      majorSigCount = majorSigCount + 1
      newSigColours[i] = baseColours[majorSigCount]

      # print(paste("searching for minors for major=", bgId, ", baseColour=", newSigColours[i], sep=''))

      # set a related colour for all other minor sigs
      minorColours = get_base_colour_extensions(newSigColours[i])
      minorCount = 1

      for(j in i+1:baSigCount)
      {
        nextSigInfo = baSigInfo[j,]
        # print(paste("testing sig: type=", nextSigInfo$Type, ", parentId=", nextSigInfo$ParentId, sep=''))

        if(!is.na(nextSigInfo$ParentId) && nextSigInfo$ParentId == bgId) # nextSigInfo$Type == "MINOR" &&
        {
          newSigColours[j] = minorColours[minorCount]
          minorCount = minorCount + 1
        }
      }
    }
  }

  # print(newSigColours)

  sigInfo = baSigInfo %>% select(Rank, BgId, Type, CancerType, Effects, SampleCount, BucketCount, RefSigs, GrpLinks, ParentId)

  sigInfo[is.na(sigInfo)] = ""
  maxTextLength = 30
  sigInfo$Effects = ifelse(stri_length(sigInfo$Effects > maxTextLength), stri_sub(sigInfo$Effects,1,maxTextLength), sigInfo$Effects)
  sigInfo$Effects = ifelse(stri_length(sigInfo$CancerType > maxTextLength), stri_sub(sigInfo$CancerType,1,maxTextLength), sigInfo$CancerType)

  sigInfo$Colour = ""
  colourCol = ncol(sigInfo)
  for(i in 1:baSigCount)
  {
    sigInfo[i,colourCol] = newSigColours[i]
  }

  sigInfo = sigInfo %>% filter(Type!="BGRD")

  evaluate_signature_fit(runType, runId, baSigs, baContribs, matrixData, summaryCounts, sampleCancerTypes,
                         bucketNames, baSigNames, plotByCancerType, viewResults, bgSigCount, printAllPlots, newSigColours, sigInfo)
}


evaluate_signature_fit<-function(runType, runId, signatures, contribution, matrixData, summaryCounts, sampleCancerTypes, bucketNames,
                            sigNamesNamed, plotByCancerType = T, viewResults = F, bgSigCount = 0, printAllPlots = T, sigColours = c(), sigInfo = data.frame())
{
  sigCount = nrow(contribution)
  hasBackgroundSigs = (bgSigCount > 0)

  print(paste("evaluating run: type=", runType, ", id=", runId, ", sigCount=", sigCount, sep=''))

  sampleNames = colnames(contribution)

  origSampleCounts = summaryCounts %>% group_by(SampleId) %>% summarise(SampleCount=sum(Count))

  sigNamesUnamed = get_signame_list(sigCount, F)
  colnames(signatures) = sigNamesUnamed

  sigNamesCombined = cbind(sigNamesUnamed, sigNamesNamed)
  colnames(sigNamesCombined) <- c("Signature", "SigName")

  print("evaluating buckets")

  bucketIndex = data.frame(as.numeric(as.character(rownames(bucketNames))))
  colnames(bucketIndex) <- c("BucketIndex")
  bucketNamesIndexed = cbind(bucketNames, bucketIndex)
  bucketNamesIndexed$BucketIndex = bucketNamesIndexed$BucketIndex-1

  # 1. Bucket Evaluation
  sigBucketData = get_bucket_data(signatures, contribution, bucketNames)
  sigBucketData = merge(sigBucketData,bucketNamesIndexed,by="Bucket",all.x=T)
  sigBucketData = merge(sigBucketData,sigNamesCombined,by="Signature",all.x=T)

  sigBucketStats = get_sig_bucket_stats(sigBucketData)

  # Signature Discovery, by looking at relative contribution of buckets
  sigBucketTopN = get_top_buckets(sigBucketData)

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
  sampleSigData = get_sig_data(signatures, contribution, sigNamesNamed, sampleNames)

  # key stats per signature
  sigStats = get_sig_stats(sampleSigData)

  # calculate and factor in residuals
  sampleSigData = append_residuals(contribution, signatures, matrixData, bucketNames, sampleSigData)
  # run again, this time bucketing samples into mutational load and cancer types

  # get cancer type and SV Count
  sampleSigData = merge(sampleSigData, sampleCancerTypes,by.x="SampleId",by.y="SampleId",all.x=TRUE)
  sampleSigData$CancerType = ifelse(is.na(sampleSigData$CancerType), 'N/A', paste(sampleSigData$CancerType, sep=""))

  # sampleSigCounts = sampleSigData %>% filter(SigName!="Residual") %>% group_by(SampleId) %>% summarise(SampleCount=sum(Count))
  sampleSigData = merge(sampleSigData, origSampleCounts, by.x="SampleId",by.y="SampleId",all.x=TRUE)

  # make the sig count total all sigs plus unalloc but not excess

  sampleSigDataNoRes = sampleSigData %>% filter(SigName!="Unalloc"&SigName!="Excess") # was 'Residuals'

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

  maxRowsPerPage = 40

  # 1. Bucket data
  title = textGrob("Bucket Summary Data & Top-N Buckets", gp=gpar(fontface="bold", fontsize=16))
  plot_bucket_summary_data(bucketSummaryData, sigBucketTopN, title, maxRowsPerPage)

  if(length(sigColours) == 0)
  {
    sigColours = get_sig_colours(sigCount)
  }

  bucketSummaryPlot = get_bucket_summary_plot(bucketSummaryData, runType)
  grid.arrange(bucketSummaryPlot, ncol = 1, nrow = 1, newpage = TRUE)

  title = textGrob("Signature-Bucket Stats", gp=gpar(fontface="bold", fontsize=16))
  pagesReqd = ceil(nrow(sigBucketStats)/maxRowsPerPage)
  for(i in 1:pagesReqd)
  {
    rowStart = (i-1)*maxRowsPerPage + 1
    rowEnd = min(rowStart+maxRowsPerPage-1,nrow(sigBucketStats))
    grid.arrange(tableGrob(sigBucketStats[rowStart:rowEnd,], rows=NULL), ncol=1, newpage = TRUE, top=title)
  }

  # 2. Default signature-bucket plot
  sigBucketsPlots = get_bucket_signatures_plot(bucketNames, signatures, sigNamesNamed)
  for(i in 1:length(sigBucketsPlots))
  {
    grid.arrange(sigBucketsPlots[[i]], ncol=1, nrow=1, newpage = TRUE)
  }

  if(printAllPlots)
  {
    # 3. Plot of Bucket Contribution for each Sample
    plot_sample_bucket_contrib(sampleBucketData, nrow(bucketNames), runType, 0) # sampleBucketTopN
  }

  # 4. Signature data
  title = textGrob("Signature Stats", gp=gpar(fontface="bold", fontsize=16))
  pagesReqd = ceil(nrow(sigStats)/maxRowsPerPage)
  for(i in 1:pagesReqd)
  {
    rowStart = (i-1)*maxRowsPerPage + 1
    rowEnd = min(rowStart+maxRowsPerPage-1,nrow(sigStats))
    grid.arrange(tableGrob(sigStats[rowStart:rowEnd,], rows = NULL), ncol=1, nrow=1, top=title, newpage = TRUE)
  }

  # add in additional colours for Excess and Residual counts
  sigColours[sigCount+1] = "black"
  sigColours[sigCount+2] = "grey30"

  # 5A. Plot sig information and a legend
  dropSigLegend = F
  if(nrow(sigInfo) > 0)
  {
    title = textGrob("Signature Info", gp=gpar(fontface="bold", fontsize=16))
    dropSigLegend = T
    # sigLegend = get_sig_plot_legend(sampleSigData, sigColours)
    grid.arrange(tableGrob(sigInfo, rows=NULL), ncol=1, nrow=1, top=title, newpage = TRUE)
    # grid.arrange(tableGrob(sigInfo, rows=NULL), sigLegend, ncol=2, nrow=1, top="Signatures", newpage = TRUE)
  }

  # 5. Relative Signature contributions by Cancer Type
  plot_cancer_sigs(sampleSigDataNoRes, sigColours, runType)

  # 6. Top 50 samples by signature, but include all other signatures as well
  plot_top_n_samples_by_sig(sampleSigDataNoRes, sigNamesNamed, 50, sigColours, runType, dropSigLegend)

  # add in additional colours for Excess and Residual counts
  # sigColours[sigCount+1] = "black"
  # sigColours[sigCount+2] = "grey30"

  # 7. Sigs with Samples by cancer type

  if(printAllPlots)
  {
    plot_sig_samples(sampleSigData, "", sigColours, runType, 0, dropSigLegend) # all samples
  }

  cancerTypes = sampleCancerTypes %>% group_by(CancerType) %>% count()

  # TEMP:
  # cancerTypes = cancerTypes %>% filter(CancerType=="Skin")

  if(plotByCancerType)
  {
    if(bgSigCount > 0)
      cancerSigColours = strip_multi_bg_colours(sigColours, bgSigCount)
    else
      cancerSigColours = sigColours

    for(cancerType in cancerTypes$CancerType)
    {
      if(!is.na(cancerType))
      {
        plot_sig_samples(sampleSigData, cancerType, cancerSigColours, runType, 0, dropSigLegend)
      }
    }
  }

  dev.off()
}




