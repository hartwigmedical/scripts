get_sig_colours<-function(sigCount = 10)
{
  if(sigCount <= 10)
  {
    sigColours = c("#ff994b","#463ec0","#d10073","#996ffb","#68b1c0","#e34bd9","#106b00","#8a392f","#98d76a","#6b3a9d","#d5c94e","#c6c0fb",
                   "#ff862c","#31528d","#d7003a","#323233","#ff4791","#01837a", "#ff748a","#777700","#ff86be")
  }
  else if(sigCount <= 20)
  {
    sigColours = c("#ff994b", "#463ec0", "#88c928", "#996ffb", "#68b1c0", "#e34bd9", "#106b00", "#d10073", "#98d76a",
                   "#6b3a9d", "#d5c94e", "#0072e2", "#ff862c", "#31528d", "#d7003a", "#323233", "#ff4791", "#01837a",
                   "#ff748a", "#777700", "#ff86be", "#4a5822", "#ffabe4", "#6a4e03", "#c6c0fb", "#ffb571", "#873659",
                   "#dea185", "#a0729d", "#8a392f")
  }
  else
  {
    sigColours = c("#ff994b", "#463ec0", "#88c928", "#996ffb", "#68b1c0", "#e34bd9", "#106b00", "#d10073", "#98d76a",
                   "#6b3a9d", "#d5c94e", "#0072e2", "#ff862c", "#31528d", "#d7003a", "#323233", "#ff4791", "#01837a",
                   "#ff748a", "#777700", "#ff86be", "#4a5822", "#ffabe4", "#6a4e03", "#c6c0fb", "#ffb571", "#873659",
                   "#dea185", "#a0729d", "#8a392f")
  }

  return (sigColours)
}

get_signame_list<-function(sigCount, asStrings)
{
  sigs = c()
  for(i in 1:sigCount)
  {
    if(asStrings & i <= 9)
    {
      sigs[i] = paste("0", i, sep='')
    }
    else
    {
      sigs[i] = i
    }
  }

  return (sigs)
}

# library(pracma)

apply_signatures<-function(matrixData, signatures, bucketCount)
{
  n_samples = dim(matrixData)[2]
  n_signatures = dim(signatures)[2]
  lsq_contribution = matrix(NA, nrow = n_signatures, ncol = n_samples)
  # lsq_reconstructed = matrix(NA, nrow = bucketCount, ncol = n_samples)

  for(i in 1:ncol(matrixData))
  {
    # print(paste("calc LSQ for i=", i, sep=''))
    y = matrixData[, i]
    y = apply(y, 1, function(x) x * 1.0) # to force it to a numeric type
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

calc_sample_residuals<-function(contribution, signatures, bucketNames, sampleBucketCounts)
{
  # contribution is sample in the columns and sig in the rows
  sigCount = ncol(signatures)
  contribTrans = t(contribution) %>% as.data.frame()
  sampleNames = colnames(contribution)
  sampleSigContribs = cbind(sampleNames, contribTrans)
  rownames(sampleSigContribs) <- NULL

  # merge with bucket-sig data
  sigContribs = cbind(signatures, bucketNames)

  # give col names so the sig coluns dont merge in the next step
  colNames = get_signame_list(sigCount,F)
  colNames[sigCount+1] = "Bucket"

  colnames(sigContribs) <- colNames

  # initially don't merge on any common fields
  samSigContribs2 = merge(sampleSigContribs, sigContribs, all.x=TRUE)
  names(samSigContribs2)[names(samSigContribs2) == 'sampleNames'] <- 'SampleId'

  # finally merge with actual bucket counts
  sampleSigContribs3 = merge(samSigContribs2, sampleBucketCounts, by.x=c("SampleId","Bucket"),by.y=c("SampleId","Bucket"),all.x=TRUE)
  sampleSigContribs3[is.na(sampleSigContribs3)] <- 0

  # convert all sig columns back to numerics
  colStart = 3
  colEnd = colStart + (sigCount*2) - 1

  # sampleSigContribs3[, c(colStart:colEnd)] <- apply(sampleSigContribs3[, c(colStart:colEnd)], 1, function(x) x * 1.0)
  sampleSigContribs3[, c(colStart:colEnd)] <- sapply(sampleSigContribs3[, c(colStart:colEnd)], as.character)
  sampleSigContribs3[, c(colStart:colEnd)] <- sapply(sampleSigContribs3[, c(colStart:colEnd)], as.numeric)

  sampleSigContribs3$SigAlloc = 0
  for(i in 1:sigCount)
  {
    sigAlloc = apply(sampleSigContribs3[,c(i+2,i+2+sigCount)], 1, function(x) x[1]*x[2])
    sampleSigContribs3$SigAlloc = sampleSigContribs3$SigAlloc + sigAlloc
  }

  sampleSigContribs3$ResidualDiff = abs(sampleSigContribs3$Count-sampleSigContribs3$SigAlloc)

  # return (sampleSigContribs3)

  sampleResidualData = (sampleSigContribs3 %>% group_by(SampleId)
                        %>% summarise(Count=sum(Count),
                                      ResidualTotal=round(sum(ResidualDiff),2),
                                      ResidualPerc=round(sum(ResidualDiff)/sum(Count),2)))

  return (sampleResidualData)
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

append_residuals<-function(contribution, signatures, bucketNames, sampleBucketCounts, sampleSigData)
{
  residuals = calc_sample_residuals(contribution, signatures, bucketNames, sampleBucketCounts)

  sampleResiduals = residuals %>% select(SampleId,ResidualTotal)
  colnames(sampleResiduals) <- c("SampleId","Count")
  sampleResiduals$Count = sampleResiduals$Count * -1.0
  sampleResiduals$SigName = "Residual"
  sampleResiduals$SigPercent = 0
  sampleResiduals$PercBucket = 0

  sampleSigData2 = rbind(sampleSigData, sampleResiduals %>% select(SampleId,SigName,SigPercent,PercBucket,Count))

  return (sampleSigData2)
}

get_sig_stats<-function(sampleSigData) {

  # key stats per signature
  sigStats = (sampleSigData %>% group_by(SigName)
              %>% summarise(SampleCount=sum(Count>0),
                            SamplePerc=round(sum(Count>0)/n_distinct(SampleId),2),
                            Count=sum(Count),
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

plot_sig_samples<-function(sampleSigData, cancerType, sigColours, varType = "SV", maxPlots = 0)
{
  cancerSigData = sampleSigData

  if(cancerType != "") {
    cancerSigData = cancerSigData %>% filter(CancerType==cancerType)
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

  if(nrow(cancerSampleSigData) > 0) {

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
    numSamples = n_distinct(cancerSampleSigData$SampleId)
    numSigs = n_distinct(cancerSampleSigData$SigName)
    samplesPerPlot = 50
    rowsPerPlot = samplesPerPlot * numSigs
    # plotCount = ceiling(numSamples/samplesPerPlot)

    rowEnd = 0
    maxRows = nrow(cancerSampleSigData)
    plotCount = 0

    while(rowEnd < maxRows)
    {
      if(plotCount == 0 & logTopNSamples)
      {
        rowStart = 1
        rowEnd = min(rowStart + topNIndex*numSigs - 1, maxRows)
      }
      else
      {
        rowStart = rowEnd + 1
        rowEnd = min(rowStart + rowsPerPlot - 1, maxRows)

        # combine stragglers onto the last graph
        if(rowEnd < maxRows & maxRows - rowEnd <= rowsPerPlot * 0.4)
        {
          # print(paste(plotCount, ": cancer=", cancerType, ", rowEnd=", rowEnd, " close to maxRows=", maxRows, sep=''))
          rowEnd = maxRows
        }
      }

      if(cancerType == "")
      {
        title = "Sig Counts by Sample"
      }
      else
      {
        title = paste("Sig Counts by Sample for ", cancerType, sep="")
      }

      sampleSigPlot <- (ggplot(cancerSampleSigData[rowStart:rowEnd,], aes(x = reorder(SampleId, -SampleCount), y = Count, fill = SigName))
      #sampleSigPlot <- (ggplot(cancerSampleSigData[rowStart:rowEnd,], aes(x = SampleId, y = Count, fill = SigName))
                        + geom_bar(stat = "identity", colour = "black")
                        + labs(x = "", y = paste(varType, " Count by Sample", sep=''))
                        + scale_fill_manual(values = sigColours)
                        + theme_bw() + theme(panel.grid.minor.x = element_blank(), panel.grid.major.x = element_blank())
                        + theme(panel.grid.minor.y = element_blank(), panel.grid.major.y = element_blank())
                        + theme(axis.text.x = element_text(angle = 90, hjust = 1,size=7)))

      if(plotCount == 0)
      {
        sampleSigPlot <- sampleSigPlot + ggtitle(title)
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

plot_top_n_samples_by_sig<-function(sampleSigData, sigNames, topN = 50, sigColours, varType = "SV")
{
  sigSamplePlots = list()
  plotIndex = 1

  # merge cancer type with sampleId
  sampleSigData = unite(sampleSigData, "SampleId", SampleId, CancerType, sep='_')

  for(sigName in sigNames)
  {
    topNSamplesBySig = head(sampleSigData %>% filter(SigName==sigName) %>% arrange(-Count),topN)
    # View(topNSamplesBySig)

    # now grab all sig data for these top-N samples
    topNSampleSigData = sampleSigData %>% filter(SampleId %in% topNSamplesBySig$SampleId)
    # View(top30SampleSigData)

    title = paste("Top Samples for Signature ", sigName, sep="")

    sampleSigPlot <- (ggplot(topNSampleSigData, aes(x = reorder(SampleId, -Count), y = Count, fill = SigName))
                      + geom_bar(stat = "identity", colour = "black")
                      + labs(x = "", y = paste(varType, " Count by Sample", sep=''))
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

calc_lsqnonneg<-function(C, d)
{
  stopifnot(is.numeric(C), is.numeric(d))
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
