# functions for handling multiple samples for the same patient
# and other combinations of variants

get_composite_id<-function(patientId,scope)
{
  compositeId = paste(patientId, ifelse(scope=="Shared","SH",ifelse(scope=="Sample1","S1","S2")), sep='_')
  return (compositeId)
}

extract_patient_id<-function(compositeId)
{
  strList = stri_split_coll(compositeId, '_')

  if(length(strList) != 1)
    return ("")

  idVec = strList[[1]]
  return (idVec[1])
}

extract_scope<-function(compositeId)
{
  strList = stri_split_coll(compositeId, '_')

  if(length(strList) != 1)
    return ("")

  idVec = strList[[1]]
  return (idVec[2])
}

# sample data is the raw sample-bucket (ie with SampleId and Bucket as fields) and must additionally include PatientId and Scope
evaluate_nmf_multiple_biopsies<-function(runType, runId, sampleData, signatures, sigNamesStr, definedBuckets, sampleCancerTypes,
                                         viewData=F, plotByCancerType=T, writeToPDF=T)
{
  print("creating compositeId")

  # 1. Create a composite ID
  sampleData$CompositeId = get_composite_id(sampleData$PatientId, sampleData$Scope)

  print("creating sample counts")

  # 2. Prepare sample-bucket counts using this composite, rather than the sample ID
  sampleCounts = sampleData %>% group_by(CompositeId,Bucket) %>% summarise(Count=sum(Count))

  # cut shared variant count in half since have been included from both samples
  sampleCounts$Count = ifelse(grepl("_SH",sampleCounts$CompositeId), sampleCounts$Count*0.5, sampleCounts$Count)

  # rename compositeId to sampleId for down-stream processing
  colnames(sampleCounts) = c("SampleId","Bucket","Count")

  print("creating matrix data")

  # 3. Create matrix data
  matrixData = sampleCounts %>% ungroup() %>% spread(SampleId, Count)
  matrixData[is.na(matrixData)] = 0

  matrixData = add_missing_buckets(matrixData, definedBuckets)

  # sort to match the defined set
  matrixData = matrixData %>% arrange(Bucket)

  bucketNames = matrixData$Bucket
  bucketNames = data.frame(bucketNames)
  colnames(bucketNames) <- c("Bucket")

  matrixData = within(matrixData, rm(Bucket))

  print("applying signatures")

  sigCount = ncol(signatures)

  # 4. Apply existing signatures
  contributions = apply_signatures(matrixData, signatures, mnvBucketCount)

  print("creating sample Id data")

  print(paste("run=", runType, ", sampleData=", nrow(sampleData), ", sigs=", sigCount, ", buckets=", nrow(bucketNames), sep=''))

  # print(paste("run=", runType, ", sampleData=", nrow(sampleData), ", sigs=", sigCount, ", buckets=", nrow(bucketNames),
  #             ", composites=", n_distinct(sampleIdData$CompositeId), ", samples=", n_distinct(sampleIdData$SampleId), sep=''))

  origSampleCounts = sampleCounts %>% group_by(SampleId) %>% summarise(OrigSampleCount=sum(Count))

  print("creating bucket data")

  bucketData = get_sample_bucket_data(matrixData, origSampleCounts, bucketNames)

  bucketData$PatientId = apply(bucketData[,c('SampleId'),drop=F], 1, function(x) extract_patient_id(x[1]))
  bucketData$Scope = apply(bucketData[,c('SampleId'),drop=F], 1, function(x) extract_scope(x[1]))

  # now work out the total variant counta across all scopes for each patient and append this to the bucket data
  patientSampleCounts = bucketData %>% group_by(PatientId) %>% summarise(PatientCount=sum(Count))
  bucketData = merge(bucketData, patientSampleCounts, by.x="PatientId", by.y="PatientId", all.x=T)

  bucketSignifDiffs = get_bucket_signif_diffs(bucketData)

  print("creating sig data")

  # repeat for signature data and plots
  sampleNames = colnames(contributions)
  sampleSigData = get_sig_data(signatures, contributions, sigNamesStr, sampleNames)
  sampleSigData = append_residuals(contributions, signatures, bucketNames, sampleCounts, sampleSigData)

  # append with patient info
  sampleSigData$PatientId = apply(sampleSigData[,c('SampleId'),drop=F], 1, function(x) extract_patient_id(x[1]))
  sampleSigData$Scope = apply(sampleSigData[,c('SampleId'),drop=F], 1, function(x) extract_scope(x[1]))

  # don't include residuals in the patient counts for signatures
  patientSampleSigCounts = sampleSigData %>% filter(SigName!="Residual") %>% group_by(PatientId) %>% summarise(PatientCount=round(sum(Count),0))

  sampleSigData = merge(sampleSigData, patientSampleSigCounts, by.x="PatientId", by.y="PatientId", all.x=T)
  sampleSigData = merge(sampleSigData, sampleCancerTypes, by.x="PatientId", by.y="PatientId", all.x=T)

  sigSignifDiffs = get_sig_signif_diffs(sampleSigData)

  if(viewData)
  {
    # View(sampleIdData)
    View(bucketData)
    View(sampleSigData)
    # View(bucketSignifDiffs)
    # View(sigSignifDiffs)
  }

  if(writeToPDF)
  {
    outputFile = paste("~/logs/r_output/pdfs/", runType, "_", runId, ".pdf", sep = "")
    print(paste("writing output to file: ", outputFile, sep=''))

    # DATA OUTPUT TO PDF
    pdf(file=outputFile, height = 14, width = 20)

    par(mar=c(1,1,1,1))

    # plot bucket data
    bucketPlots = get_sample_bucket_mult_biop_plots(bucketData, nrow(bucketNames), runType, 0)

    plot_matching_graphs(bucketPlots)

    if(nrow(bucketSignifDiffs)>0)
    {
      title = textGrob("Buckets with Significant Private Increase", gp=gpar(fontface="bold", fontsize=16))
      grid.arrange(tableGrob(head(bucketSignifDiffs,30), rows = NULL), ncol = 1, nrow = 1, top=title, newpage = TRUE)
    }

    # plot sig data
    sigColours = get_sig_colours(sigCount)
    sigColoursPlusRes = sigColours
    sigColoursPlusRes[sigCount+1] = "#000000"

    sigPlots = get_sig_samples_mult_biop_plots(sampleSigData, "", sigColoursPlusRes, runType, 0) # all samples
    plot_matching_graphs(sigPlots)

    if(nrow(sigSignifDiffs)>0)
    {
      title = textGrob("Signatures with Significant Private Increase", gp=gpar(fontface="bold", fontsize=16))
      grid.arrange(tableGrob(head(sigSignifDiffs,30), rows = NULL), ncol = 1, nrow = 1, top=title, newpage = TRUE)
    }

    cancerTypes = sampleCancerTypes %>% group_by(CancerType) %>% count()

    if(plotByCancerType)
    {
      for(cancerType in cancerTypes$CancerType)
      {
        sigPlots = get_sig_samples_mult_biop_plots(sampleSigData, cancerType, sigColoursPlusRes, runType, 0)

        if(length(sigPlots) > 0)
          plot_matching_graphs(sigPlots)
      }
    }

    # plot key differences

    dev.off()
  }
}


get_sample_bucket_mult_biop_plots<-function(sampleBucketData, bucketCount, varType = "SV", maxPlots = 4)
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

  if(nrow(sbData) > 0)
  {
    sbData = sbData %>% arrange(-PatientCount,PatientId,Bucket)

    sbDataSH = sbData %>% filter(Scope=="SH")
    sbDataS1 = sbData %>% filter(Scope=="S1")
    sbDataS2 = sbData %>% filter(Scope=="S2")

    bucketPlots = list()

    if(nrow(sbDataSH) != nrow(sbDataS1) | nrow(sbDataS1) != nrow(sbDataS2))
    {
      print("WARN: unequal bucket data set counts")
      return (bucketPlots)
    }

    plotIndex = 1

    samplesPerPlot = 100
    rowsPerPlot = samplesPerPlot * bucketCount

    rowEnd = 0
    maxRows = nrow(sbDataSH)
    plotCount = 0

    while(rowEnd < maxRows)
    {
      if(plotCount == 0)
      {
        rowStart = 1
        rowEnd = min(rowStart + (samplesPerPlot*0.75)*bucketCount - 1, maxRows) # less to leave space for the legend
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

      title = "Bucket % by Sample"

      for(plotIndex in 1:3)
      {
        # print(paste("plotCount=", plotCount, ", index=", plotIndex, ", rowStart=", rowStart, ", rowEnd=", rowEnd, sep=''))
        includeInfo = F

        if(plotIndex == 1)
        {
          includeInfo = T
          sbPlotData = sbDataSH
        }
        else if(plotIndex == 2)
        {
          sbPlotData = sbDataS1
        }
        else
        {
          sbPlotData = sbDataS2
        }

        plot = (ggplot(sbPlotData[rowStart:rowEnd,], aes(x = reorder(PatientId, -PatientCount), y=Percent, fill=Bucket))
                + geom_bar(stat = "identity", colour = "black", size=0.01)
                + labs(x ="", y="Bucket Percentages")
                + scale_fill_manual(values = bucketColours)
                + theme_bw() + theme(panel.grid.minor.x = element_blank(), panel.grid.major.x = element_blank())
                + theme(panel.grid.minor.y = element_blank(), panel.grid.major.y = element_blank())
                + theme(axis.text.x = element_text(angle = 90, hjust = 1,size=7)))

        if(includeInfo)
        {
          plot <- plot + ggtitle(title) + theme(legend.text = element_text(size=7))
        }
        else
        {
          # remove patient IDs after the first plot
          plot <- plot + theme(axis.text.x=element_blank())
        }

        if(plotCount > 0)
        {
          # remove legend after the first plot
          plot <- plot + theme(legend.position="none")
        }

        plotCount = plotCount + 1
        bucketPlots[[plotCount]] = plot
      }

      if(maxPlots > 0 && plotCount >= maxPlots)
        break
    }
  }

  return (bucketPlots)
}

get_sig_samples_mult_biop_plots<-function(sampleSigData, cancerType, sigColours, varType = "SV", maxPlots = 0)
{
  if(cancerType != "")
  {
    cancerSSData = sampleSigData %>% filter(CancerType==cancerType)
  }
  else
  {
    cancerSSData = sampleSigData
  }

  sigPlots = list()
  plotIndex = 1

  if(nrow(cancerSSData) == 0)
  {
    return (sigPlots)
  }

  cancerSSData = cancerSSData %>% arrange(-PatientCount,PatientId,SigName)

  cssDataSH = cancerSSData %>% filter(Scope=="SH")
  cssDataS1 = cancerSSData %>% filter(Scope=="S1")
  cssDataS2 = cancerSSData %>% filter(Scope=="S2")

  if(nrow(cssDataSH) != nrow(cssDataS1) | nrow(cssDataS1) != nrow(cssDataS2))
  {
    print(paste("WARN: unequal sig data set counts for cancerType=", cancerType, sep=''))
    return (sigPlots)
  }

  samplesPerPlot = 100
  sigCount = n_distinct(cancerSSData$SigName)
  rowsPerPlot = samplesPerPlot * sigCount

  rowEnd = 0
  maxRows = nrow(cssDataSH)
  plotCount = 0

  # plot the first X samples separately if they differ substantiately
  topNLevel = 6
  topNIndex = topNLevel * sigCount

  logTopNSamples = F
  if(nrow(cssDataSH) > topNIndex)
  {
    maxCount = max(cssDataSH$PatientCount)
    nthCount = min(head(cssDataSH,topNIndex)$PatientCount)

    if(maxCount > 3 * nthCount)
    {
      logTopNSamples = T
      print(paste("splitting top samples: maxCount=", maxCount, ", nthCount=", nthCount, sep=''))
    }
  }

  # print(paste("sigCount=", sigCount, ", rowsPerPlot=", rowsPerPlot, ", maxRows=", maxRows, sep=''))

  while(rowEnd < maxRows)
  {
    if(plotCount == 0)
    {
      rowStart = 1

      if(logTopNSamples)
      {
        rowEnd = min(rowStart + (topNLevel-1)*sigCount - 1, maxRows)
      }
      else
      {
        rowEnd = min(rowStart + (samplesPerPlot*0.75)*sigCount - 1, maxRows) # less to leave space for the legend
      }
    }
    else
    {
      rowStart = rowEnd + 1
      rowEnd = min(rowStart + rowsPerPlot - 1, maxRows)

      # combine stragglers onto the last graph
      if(rowEnd < maxRows & maxRows - rowEnd < rowsPerPlot * 0.4)
      {
        rowEnd = maxRows
      }
    }

    if(cancerType == "")
      title = "Sig Counts by Sample"
    else
      title = paste("Sig Counts by Sample for ", cancerType, sep="")

    for(plotIndex in 1:3)
    {
      # print(paste("plotCount=", plotCount, ", index=", plotIndex, ", rowStart=", rowStart, ", rowEnd=", rowEnd, sep=''))
      includeInfo = F

      if(plotIndex == 1)
      {
        includeInfo = T
        cssPlotData = cssDataSH
        plotTitle = "Shared"
      }
      else if(plotIndex == 2)
      {
        cssPlotData = cssDataS1
        plotTitle = "Sample 1"
      }
      else
      {
        cssPlotData = cssDataS2
        plotTitle = "Sample 2"
      }

      plot = (ggplot(cssPlotData[rowStart:rowEnd,], aes(x = reorder(PatientId, -PatientCount), y=Count, fill=SigName))
              + geom_bar(stat = "identity", colour = "black", size=0.01)
              + labs(x ="", y="Variant Counts")
              + scale_fill_manual(values = sigColours)
              + theme_bw() + theme(panel.grid.minor.x = element_blank(), panel.grid.major.x = element_blank())
              + theme(panel.grid.minor.y = element_blank(), panel.grid.major.y = element_blank())
              + theme(axis.text.x = element_text(angle = 90, hjust = 1,size=7)))

      if(includeInfo)
      {
        plot <- plot + ggtitle(title) + theme(legend.text = element_text(size=7))
      }
      else
      {
        # remove patient IDs after the first plot
        plot <- plot + theme(axis.text.x=element_blank())
      }

      if(plotCount > 0)
      {
        # remove legend after the first plot
        plot <- plot + ggtitle(plotTitle) + theme(legend.position="none")
      }

      plotCount = plotCount + 1
      sigPlots[[plotCount]] = plot
    }

    if(maxPlots > 0 && plotCount >= maxPlots)
      break
  }

  return (sigPlots)
}

plot_matching_graphs<-function(plotList, firstHasLegend=T)
{
  plotCount = length(plotList)

  if(plotCount == 0)
  {
    print("No plots provided")
  }

  i = 1
  pageCount = 1
  while(i < plotCount)
  {
    # print(paste("Plotting page ", pageCount, ", index=", i, " of ", plotCount, sep=''))

    if(pageCount == 1 & firstHasLegend)
    {
      plot1 = plotList[[1]]
      plotLegend <- get_legend(plot1)
      plot1 = plot1 + theme(legend.position="none")

      grid.arrange(plot_grid(plot1, plotList[[i+1]], plotList[[i+2]], ncol=1),
                   plotLegend,
                   ncol=2, widths = c(0.8, 0.2), newpage=TRUE)

      # only include the legend on the first page (in a separate column)
      # ggdraw(plot_grid(plot_grid(plot1, plotList[[i+1]], plotList[[i+2]], ncol=1, align='v'),
      #                  plot_grid(plotLegend, ncol=1),
      #                  rel_widths=c(0.8, 0.2)))
    }
    else
    {
      grid.arrange(plot_grid(plotList[[i]], plotList[[i+1]], plotList[[i+2]], ncol=1), newpage=T)
      #ggdraw(plot_grid(plotList[[i]], plotList[[i+1]], plotList[[i+2]], ncol=1))
    }

    pageCount = pageCount + 1
    i = i+3
  }
}


get_bucket_signif_diffs<-function(bucketData, sharedPercentMin = 0.05, sigDiffPerc = 0.1)
{
  bucketPerPatient = (bucketData %>% group_by(PatientId,Bucket)
                      %>% summarise(TotalCount=sum(Count),
                                    Shared_Perc=round(sum(ifelse(Scope=="SH",Percent,0)),2),
                                    S1_Perc=round(sum(ifelse(Scope=="S1",Percent,0)),2),
                                    S2_Perc=round(sum(ifelse(Scope=="S2",Percent,0)),2),
                                    Max_Perc=round(max(Percent),2)))

  bucketPerPatient$MaxvsSharedPerc = bucketPerPatient$Max_Perc-bucketPerPatient$Shared_Perc

  bucketSignifDiffs = (bucketPerPatient
                       %>% filter((S1_Perc-Shared_Perc>=sigDiffPerc|S2_Perc-Shared_Perc>=sigDiffPerc)&Shared_Perc>=sharedPercentMin)
                       %>% arrange(-MaxvsSharedPerc))

  return (bucketSignifDiffs)
}

get_sig_signif_diffs<-function(sigData, sharedPercentMin = 0.05, sigDiffPerc = 0.1)
{
  sigsPerPatient = (sigData %>% filter(SigName!="Residual") %>% group_by(PatientId,SigName)
                    %>% summarise(CancerType=first(CancerType),
                                  TotalCount=round(sum(Count),0),
                                  Shared_Perc=round(sum(ifelse(Scope=="SH",SigPercent,0)),2),
                                  S1_Perc=round(sum(ifelse(Scope=="S1",SigPercent,0)),2),
                                  S2_Perc=round(sum(ifelse(Scope=="S2",SigPercent,0)),2),
                                  Max_Perc=round(max(SigPercent),2)))

  sigsPerPatient$MaxvsSharedPerc = sigsPerPatient$Max_Perc-sigsPerPatient$Shared_Perc

  sigSignifDiffs = (sigsPerPatient
                    %>% filter((S1_Perc-Shared_Perc>=sigDiffPerc|S2_Perc-Shared_Perc>=sigDiffPerc)&Shared_Perc>=sharedPercentMin)
                    %>% arrange(-MaxvsSharedPerc))

  return (sigSignifDiffs)
}

