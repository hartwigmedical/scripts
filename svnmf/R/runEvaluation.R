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

trim_ba_sig_names<-function(baSigNames)
{
  baSigCount = length(baSigNames)

  baSigStrList = get_signame_list(baSigCount, T)
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

  return (baSigNames)
}

produce_signature_report<-function(runType, runId, sigFile, contribFile, sigInfoFile, sigAllocFile,
                                   matrixData, summaryCounts, sampleCancerTypes, bucketNames,
                                   plotByCancerType = T, viewResults = F, printAllPlots = T)
{
  baContribs = as.matrix(read.csv(file=contribFile, stringsAsFactors=F))
  baSigs = as.matrix(read.csv(file=sigFile, stringsAsFactors=F))
  sigAllocs = as.data.frame(read.csv(file=sigAllocFile, stringsAsFactors=F))
  colnames(sigAllocs) = c("Signature", "BgId", "SampleId", "Count", "SigPercent")

  # TEMP until next java run
  sigAllocs$Signature = ifelse(sigAllocs$BgId=='Excess',42,sigAllocs$Signature)
  sigAllocs$Count = ifelse(sigAllocs$BgId=='Excess',-sigAllocs$Count,sigAllocs$Count) # negate the excess counts

  baSigNames = colnames(baSigs)
  baSigCount = ncol(baSigs)

  baSigNumList = get_signame_list(baSigCount, F)
  colnames(baSigs) <- baSigNumList

  # trim sig names
  baSigNames = trim_ba_sig_names(baSigNames)

  baSigInfo = as.data.frame(read.csv(file=sigInfoFile, stringsAsFactors=F))

  # look for links between the sigs
  baSigCount = nrow(baSigInfo)

  baseColours = get_base_colours()
  generalColours = get_unique_colours()

  # look for links between the sigs
  bgSigCount = 0
  newSigColours = c()

  # set all background sigs to grey
  for(i in 1:baSigCount)
  {
    sigInfo = baSigInfo[i,]

    if(sigInfo$Type == "BGRD")
    {
      bgSigCount = bgSigCount + 1
      newSigColours[i] = "grey60"
    }
    else
    {
      newSigColours[i] = "unset" # unset for now
    }
  }

  # set specific colours for majors and minors
  majorSigCount = 0
  for(i in 1:baSigCount)
  {
    sigInfo = baSigInfo[i,]

    if(sigInfo$Type == "MAJOR" && (is.na(sigInfo$ParentId) || sigInfo$ParentId==""))
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

  # set colours for other misc sigs
  otherSigCount = 0
  for(i in 1:baSigCount)
  {
    if(newSigColours[i] == "unset")
    {
      otherSigCount = otherSigCount + 1
      newSigColours[i] = generalColours[otherSigCount]
    }
  }

  sigInfo = baSigInfo %>% select(Rank, BgId, Type, Discovery, CancerType, Effects, SampleCount, BucketCount, RefSigs, GrpLinks, ParentId)

  sigInfo[is.na(sigInfo)] = ""
  maxTextLength = 30
  sigInfo$Effects = ifelse(stri_length(sigInfo$Effects > maxTextLength), stri_sub(sigInfo$Effects,1,maxTextLength), sigInfo$Effects)
  sigInfo$CancerType = ifelse(stri_length(sigInfo$CancerType > maxTextLength), stri_sub(sigInfo$CancerType,1,maxTextLength), sigInfo$CancerType)

  sigInfo$Colour = ""
  colourCol = ncol(sigInfo)
  for(i in 1:baSigCount)
  {
    sigInfo[i,colourCol] = newSigColours[i]
  }

  sigInfo = sigInfo %>% filter(Type!="BGRD")

  evaluate_signature_fit(runType, runId, baSigs, baContribs, matrixData, summaryCounts, sampleCancerTypes,
                         bucketNames, baSigNames, plotByCancerType, viewResults, bgSigCount, printAllPlots, newSigColours, sigInfo, sigAllocs)
}


evaluate_signature_fit<-function(runType, runId, signatures, contribution, matrixData, summaryCounts, sampleCancerTypes, bucketNames,
                            sigNamesNamed, plotByCancerType = T, viewResults = F, bgSigCount = 0, printAllPlots = T, sigColours = c(),
                            sigInfo = data.frame(), sigAllocs = data.frame())
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

  if(nrow(sigAllocs) > 0)
  {
    sigNamesCombined = rbind(sigNamesCombined, c(sigCount+1,"Unalloc"))
    sigNamesCombined = rbind(sigNamesCombined, c(sigCount+2,"Excess"))

    # ensure an entry in every sig or every sample
    sampleSigFullSet = merge(sampleNames, sigNamesCombined)
    colnames(sampleSigFullSet) = c("SampleId", "Signature", "SigName")

    sampleSigData = merge(sigAllocs, sampleSigFullSet, by=c("SampleId", "Signature"), all=T)
    sampleSigData[is.na(sampleSigData)] = 0
    sampleSigData$PercBucket = round(sampleSigData$SigPercent/0.1)*0.1
    sampleSigData = within(sampleSigData, rm(BgId))
  }
  else
  {
    sampleSigData = get_sig_data(signatures, contribution, sigNamesNamed, sampleNames)

    # calculate and factor in residuals
    sampleSigData = append_residuals(contribution, signatures, matrixData, bucketNames, sampleSigData)
  }

  # get cancer type and SV Count
  sampleSigData = merge(sampleSigData, sampleCancerTypes,by.x="SampleId",by.y="SampleId",all.x=TRUE)
  sampleSigData$CancerType = ifelse(is.na(sampleSigData$CancerType), 'N/A', paste(sampleSigData$CancerType, sep=""))
  sampleSigData = merge(sampleSigData, origSampleCounts, by.x="SampleId",by.y="SampleId",all.x=TRUE)
  sampleSigDataNoResiduals = sampleSigData %>% filter(SigName!="Unalloc"&SigName!="Excess")

  # key stats per signature
  sigStats = get_sig_stats(sampleSigDataNoResiduals)

  # make the sig count total all sigs plus unalloc but not excess
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

  if(bgSigCount > 0)
    cancerSigColours = strip_multi_bg_colours(sigColours, bgSigCount)
  else
    cancerSigColours = sigColours

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
  plot_cancer_sigs(sampleSigDataNoResiduals, sigColours, runType)

  # 6. Top 50 samples by signature, but include all other signatures as well
  plot_top_n_samples_by_sig(sampleSigDataNoResiduals, sigNamesNamed, 50, sigColours, runType, dropSigLegend)

  # 7. Sigs with Samples by cancer type

  if(printAllPlots)
  {
    plot_sig_samples(sampleSigData, "", sigColours, runType, 0, dropSigLegend) # all samples
  }

  cancerTypes = sampleCancerTypes %>% group_by(CancerType) %>% count()

  # TEMP:
  cancerTypes = cancerTypes %>% filter(CancerType=="Skin")

  if(plotByCancerType)
  {
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




