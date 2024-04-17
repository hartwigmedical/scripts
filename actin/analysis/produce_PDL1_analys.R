library(tidyr)
library(dplyr)
library(tidyverse)
library(ggplot2)
library(grid)
library(gridExtra)
library(DBI)
library(RMySQL)

wd <- paste0(Sys.getenv("HOME"), "/hmf/tmp/")
dbConnect <- dbConnect(MySQL(), dbname='hmfpatients', groups="RAnalysis")

# CD274 DNA amp count vs total number of samples
ampCountQuery <- "
select count(dr.driver = 'AMP') as ampCount, totalCount
from driverCatalog dr
inner join datarequest da on da.sampleId = dr.sampleId 
inner join hpc on hpc.sampleId = dr.sampleId
inner join (select count(da.sampleId) as totalCount from datarequest da inner join hpc on da.sampleId=hpc.sampleId) total
where dr.gene = 'CD274';
"
ampCountResult <- dbGetQuery(dbConnect, ampCountQuery)
ampPercentage <- round(ampCountResult$ampCount/ampCountResult$totalCount*100,2)
print(paste("CD274amp is found in", ampCountResult$ampCount, "samples in the Hartwig database (total", ampCountResult$totalCount, ",", ampPercentage,"%)"))

## Number of ACTIN samples (full cohort) with CD274 amp
ACTNampCountQuery <- "
select count(driver = 'AMP') as ampCount
from driverCatalog 
where sampleId like 'ACTN%' and gene = 'CD274';
"
ACTNampCountResult <- dbGetQuery(dbConnect, ACTNampCountQuery)
print(paste("Number of ACTIN samples with CD274 amplification:", ACTNampCountResult))

## TPM distribution in ACTIN samples
tpmQuery <- "
select tpm from geneExpression e
inner join datarequest d on d.sampleId = e.sampleId
where e.sampleId like 'ACTN%' 
and gene = 'CD274';
"
tmpResult <- dbGetQuery(dbConnect, tpmQuery)
tmpValues <- unlist(tmpResult)
tmp <- as.numeric(tmpValues)

percentile_90 <- quantile(tmp, 0.9)

hist(tmp, breaks = seq(0, max(tmp) + 1, by = 1), freq = TRUE, col = "lightblue", main = "CD274 transcripts per million in ACTIN samples",
     xlab = "TPM", ylab = "Number of Samples")
abline(v = percentile_90, col = "red", lwd = 2)

## top 10 percent CD274 tpm in entire HMF db
tpmTop10Query <- "
select d.sampleId
from geneExpression g
inner join datarequest d on d.sampleId = g.sampleId
inner join hpc on hpc.sampleId = g.sampleId
inner join rnaStatistics r on g.sampleId=r.sampleId
where gene = 'CD274' and r.qcStatus='PASS'
and g.percentileCohort >= 0.9;
"
tpmTop10Result <- dbGetQuery(dbConnect, tpmTop10Query)

## CD274 AMP on DNA level in top 10% CD274 expression samples
CD274ampQuery <- "
select dr.sampleId
from driverCatalog dr
inner join datarequest da on da.sampleId = dr.sampleId
inner join hpc on hpc.sampleId = dr.sampleId
inner join rnaStatistics r on dr.sampleId=r.sampleId
where dr.gene = 'CD274' and r.qcStatus='PASS' and dr.driver = 'AMP';
"
CD274ampResult <- dbGetQuery(dbConnect, CD274ampQuery)

CD274ampInTopTPMResult <- tpmTop10Result %>% inner_join(CD274ampResult, by = "sampleId") %>% nrow()

ampAndExpressionQuery <- "
select count(dr.driver = 'AMP') as ampCount
from driverCatalog dr
inner join datarequest da on da.sampleId = dr.sampleId
inner join hpc on hpc.sampleId = dr.sampleId
inner join geneExpression g on g.sampleID = dr.sampleId
inner join rnaStatistics r on dr.sampleId=r.sampleId
where dr.gene = 'CD274' and r.qcStatus='PASS' and g.gene = 'CD274';
"
ampAndExpressionResult <- dbGetQuery(dbConnect, ampAndExpressionQuery)

## Percentage of samples with CD274amp on DNA level in top 10% CD274 expression samples
CD274ampPercIntop10TPM <- CD274ampInTopTPMResult/count(tpmTop10Result)*100
print(paste("Of top 10% CD274 TPM samples", "(total:", count(tpmTop10Result), ")" , "CD274amp is found in:", CD274ampInTopTPMResult, ",", CD274ampPercIntop10TPM, "%"))
print(paste("Total number of CD274 amplification samples (with RNA Expression data) is:", ampAndExpressionResult))
percCD274ampWithExpressionAbove90thPercentile <- CD274ampInTopTPMResult/ampAndExpressionResult*100
print(paste(percCD274ampWithExpressionAbove90thPercentile, "% of CD274 amp samples are above the 90th percentile of CD274 Gene Expression (RNA)", "(", CD274ampInTopTPMResult, "of", ampAndExpressionResult, ")"))

## Copy number (DNA) and percentileCohort (RNA expr) plot for CD274
copyNumAndExprPercQuery <- "
select ex.sampleId, minCopyNumber, percentileCohort 
from geneCopyNumber cn
inner join datarequest dr on dr.sampleId = cn.sampleId
inner join geneExpression ex on ex.sampleId = cn.sampleId
inner join hpc on hpc.sampleId = dr.sampleId
inner join rnaStatistics r on dr.sampleId=r.sampleId
where ex.gene = 'CD274' and r.qcStatus='PASS' and cn.gene = 'CD274';
"
copyNumAndExprPercResult <- dbGetQuery(dbConnect, copyNumAndExprPercQuery)

##CD274 AMP samples with RNA data
AMPwithRNAQuery <- "
select dr.sampleId 
from driverCatalog dr
inner join datarequest da on da.sampleId = dr.sampleId
inner join geneExpression g on g.sampleId = dr.sampleId
inner join hpc on hpc.sampleId = dr.sampleId
inner join rnaStatistics r on dr.sampleId=r.sampleId
where dr.gene = 'CD274' and r.qcStatus='PASS' and g.gene = 'CD274' and driver = 'AMP';
"
AMPwithRNAResult <- dbGetQuery(dbConnect, AMPwithRNAQuery)

## Plot of Copy number (DNA) versus percentile cohort for CD274 expression - DNA amp samples in red
ggplot(copyNumAndExprPercResult, aes(x = percentileCohort, y = minCopyNumber)) +
  geom_point(data = subset(copyNumAndExprPercResult, sampleId %in% AMPwithRNAResult$sampleId), 
             aes(color = "Samples with amplification on DNA level"), 
             size = 3) +
  geom_point(data = subset(copyNumAndExprPercResult, !(sampleId %in% AMPwithRNAResult$sampleId)), 
             aes(color = "No amplification"), 
             size = 1) +
  labs(x = "RNA Expression - Percentile Cohort", y = "DNA - (Minimum) Copy Number") +
  ggtitle("CD274: Copy Number vs. RNA percentile cohort") +
  theme_minimal() +
  scale_color_manual(values = c("No amplification" = "black", "Samples with amplification on DNA level" = "red"),
                     labels = c( "No amplification", "Samples with amplification on DNA level"),
                     name = "Legend")


## RNA expression of CD274 in samples with MYC amp and NO CD274 amp
MYCampQuery <- "
select d.sampleId, d.driver, d.gene as DNAgene, g.gene as RNAgene, g.tpm, g.percentileCohort, cn.minCopyNumber 
from driverCatalog d
inner join datarequest da on da.sampleId = d.sampleId
inner join hpc on hpc.sampleId = d.sampleId
inner join geneExpression g on g.sampleId = d.sampleId and g.gene = 'CD274'
inner join geneCopyNumber cn on cn.sampleId = d.sampleId and cn.gene = 'CD274'
inner join rnaStatistics r on da.sampleId=r.sampleId
where r.qcStatus='PASS' and d.gene = 'MYC'
and d.sampleId not in (
  select sampleId from driverCatalog where gene = 'CD274' and d.driver = 'AMP'
)
and d.driver = 'AMP'
group by d.sampleId;
"

MYCampQueryResult <- dbGetQuery(dbConnect, MYCampQuery)
MYCampCD274ExprAbovePercentile90 <- tpmTop10Result %>% inner_join(MYCampQueryResult, by = "sampleId") %>% nrow()
print(paste(MYCampCD274ExprAbovePercentile90, "of", count(MYCampQueryResult), "(", MYCampCD274ExprAbovePercentile90/count(MYCampQueryResult) *100, "%) MYC amp samples are above the 90th percentile of CD274 expression"))

## Plot of MYC Copy number (DNA) versus percentile cohort for CD274 expression 
ggplot(MYCampQueryResult, aes(x = percentileCohort, y = minCopyNumber)) +
  geom_point(size = 3) +
  labs(x = "CD274 RNA Expression - Percentile Cohort", y = "MYC DNA - (Minimum) Copy Number") +
  ggtitle("MYC Min Copy Number vs. CD274 RNA percentile cohort") +
  theme_minimal()

## Copy number (DNA) and percentileCohort (RNA expr) plot for JAK2
JAK2copyNumAndExprPercQuery <- "
select d.sampleId, d.driver, d.gene as DNAgene, g.gene as RNAgene, g.tpm, g.percentileCohort, cn.minCopyNumber 
from driverCatalog d
inner join datarequest da on da.sampleId = d.sampleId
inner join hpc on hpc.sampleId = d.sampleId
inner join geneExpression g on g.sampleId = d.sampleId and g.gene = 'CD274'
inner join geneCopyNumber cn on cn.sampleId = d.sampleId and cn.gene = 'JAK2'
inner join rnaStatistics r on da.sampleId=r.sampleId
where r.qcStatus='PASS' and d.gene = 'JAK2'
and d.sampleId not in (
    select sampleId from driverCatalog where gene = 'CD274' and d.driver = 'AMP'
)
and d.driver = 'AMP'
group by d.sampleId;
"
JAK2copyNumAndExprPercResult <- dbGetQuery(dbConnect, JAK2copyNumAndExprPercQuery)
JAK2ampInTop10Expr <- tpmTop10Result %>% inner_join(JAK2copyNumAndExprPercResult, by = "sampleId") %>% nrow()
print(paste(JAK2ampInTop10Expr, "of", count(JAK2copyNumAndExprPercResult), "samples with JAK2 amp (and no CD274 amp) have a CD274 expression above the 90th percentile"))

## Plot of JAK2 Copy number (DNA) versus percentile cohort for CD274 expression 
ggplot(JAK2copyNumAndExprPercResult, aes(x = percentileCohort, y = minCopyNumber)) +
  geom_point(size = 3) +
  labs(x = "CD274 RNA Expression - Percentile Cohort", y = "JAK2 DNA - (Minimum) Copy Number") +
  ggtitle("JAK2 Copy Number vs. CD274 RNA percentile cohort") +
  theme_minimal()

## RNA expression in samples with CD274 in-frame fusion
CD274FusionQuery <- "
select g.sampleId, gene, tpm, percentileCancer, percentileCohort from geneExpression g
inner join datarequest d on d.sampleId = g.sampleId
inner join hpc on hpc.sampleId = g.sampleId
inner join rnaStatistics r on d.sampleId=r.sampleId
where g.sampleId in (select sampleId from svFusion where name like '%_CD274%' and phased = 'INFRAME')
and r.qcStatus='PASS'
and g.gene = 'CD274';
"
CD274FusionResult <- dbGetQuery(dbConnect, CD274FusionQuery)

CD274FusionInTop10Expr <- tpmTop10Result %>% inner_join(CD274FusionResult, by = "sampleId") %>% nrow()
percentage_HighExpr <- CD274FusionInTop10Expr / count(tpmTop10Result) * 100
print(paste(CD274FusionInTop10Expr, "of", count(CD274FusionResult), "samples with CD274 Fusion have a CD274 expression above the 90th percentile (", CD274FusionInTop10Expr/count(CD274FusionResult) *100, "% )"))
print(paste("Of top 10% CD274 TPM samples", "(total:", count(tpmTop10Result), ")" , "a CD274 Fusion is found in:", CD274FusionInTop10Expr, ",", round(percentage_HighExpr, 2), "%"))

## RNA expression in samples with CD274 3UTR disruption
CD274disruption3UTRQuery <- "
select g.sampleId, gene, tpm, percentileCancer, percentileCohort from geneExpression g
inner join datarequest d on d.sampleId = g.sampleId
inner join hpc on hpc.sampleId = g.sampleId
inner join rnaStatistics r on d.sampleId=r.sampleId
where g.sampleId in (select sampleId from svBreakend where gene = 'CD274' and codingContext = 'UTR_3P')
and r.qcStatus='PASS'
and gene = 'CD274';
"
CD274disruption3UTRResult <- dbGetQuery(dbConnect, CD274disruption3UTRQuery)
CD274disruption3UTRInTop10Expr <- tpmTop10Result %>% inner_join(CD274disruption3UTRResult, by = "sampleId") %>% nrow()
print(paste(CD274disruption3UTRInTop10Expr, "of", count(CD274disruption3UTRResult), "samples with CD274 3UTR disruption have a CD274 expression above the 90th percentile (", CD274disruption3UTRInTop10Expr/count(CD274disruption3UTRResult) *100, "% )"))

## RNA expression in samples with CD274 disruptive structural variant
CD274disruptiveSVQuery <- "
select g.sampleId, g.tpm, g.percentileCohort 
from geneExpression g 
inner join datarequest d on d.sampleId=g.sampleId
inner join hpc on hpc.sampleId = g.sampleId
inner join rnaStatistics r on d.sampleId=r.sampleId
where g.sampleId in (select sampleId from svBreakend where gene='CD274' and disruptive and codingContext != 'UTR_3P')
and r.qcStatus='PASS'
and gene='CD274' order by g.percentileCohort desc;
"
CD274disruptiveSVResult <- dbGetQuery(dbConnect, CD274disruptiveSVQuery)
CD274disruptiveSVResultInTop10Expr <- tpmTop10Result %>% inner_join(CD274disruptiveSVResult, by = "sampleId") %>% nrow()
print(paste(CD274disruptiveSVResultInTop10Expr, "of", count(CD274disruptiveSVResult), "samples with CD274 disruptive structural variant have a CD274 expression above the 90th percentile (", CD274disruptiveSVResultInTop10Expr/count(CD274disruptiveSVResult) *100, "% )"))

## RNA expression in samples with either CD274 amp or CD274 inframe fusion
CD274ampOrFusionQuery <- "
select da.sampleId, g.gene as RNAgene, g.tpm, g.percentileCohort 
from datarequest da 
inner join hpc on hpc.sampleId = da.sampleId
inner join geneExpression g on g.sampleId = da.sampleId and g.gene = 'CD274'
inner join rnaStatistics r on da.sampleId=r.sampleId
where (
  da.sampleId in (select sampleId from driverCatalog where gene='CD274' and driver='AMP')
  or da.sampleId in (select sampleId from svFusion where name like '%_CD274%' and phased = 'INFRAME')
)
and r.qcStatus='PASS';
"
CD274ampOrFusionResult <- dbGetQuery(dbConnect, CD274ampOrFusionQuery)
CD274ampOrFusionResultInTop10Expr <- tpmTop10Result %>% inner_join(CD274ampOrFusionResult, by = "sampleId") %>% nrow()
print(paste(CD274ampOrFusionResultInTop10Expr, "of", count(CD274ampOrFusionResult), "samples with CD274 amp or fusion have a CD274 expression above the 90th percentile (", CD274ampOrFusionResultInTop10Expr/count(CD274ampOrFusionResult) *100, "% )"))

CD274ampOrFusionResultAsPartOfTop10Expr <- CD274ampOrFusionResultInTop10Expr / count(tpmTop10Result) * 100
print(paste(CD274ampOrFusionResultAsPartOfTop10Expr, "% of samples above 90th percentile expression have either CD274 amp or fusion"))

## Combine CD274amp + fusion + MYCamp + 3UTR
combinedResults <- CD274ampOrFusionResult %>% 
  full_join(CD274disruptiveSVResult, by = "sampleId") %>%
  full_join(CD274disruption3UTRResult, by = "sampleId")

length(combinedResults$sampleId)
combinedResultsInTop10Expr <- tpmTop10Result %>% inner_join(combinedResults, by = "sampleId") %>% nrow()
combinedInTopPercentage <- combinedResultsInTop10Expr / count(combinedResults) * 100
combinedInTopPercentage
combinedAsPartOfTopPercentage <- combinedResultsInTop10Expr / count(tpmTop10Result) * 100
combinedAsPartOfTopPercentage

print(paste(combinedResultsInTop10Expr, "of", length(combinedResults$sampleId), "samples with CD274 amp, fusion, 3UTR disruption or other disruptive SV have a CD274 expression above the 90th percentile (", combinedResultsInTop10Expr / length(combinedResults$sampleId) * 100, "% )"))
print(paste(combinedResultsInTop10Expr, "of", length(tpmTop10Result$sampleId), "(", combinedAsPartOfTopPercentage, "%)", "of samples above 90th percentile expression have either CD274 amp, fusion, 3UTR disruption or other disruptive SV"))

dbDisconnect(dbConnect)

