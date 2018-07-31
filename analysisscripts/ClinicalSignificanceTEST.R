library(dplyr)
library(tidyr)
library(RMySQL)
library(data.table)

get_rano_data<-function(dbConnect,sampleString="")
{
  query = paste("SELECT * FROM ranoMeasurement ",sep="")
  
  if (sampleString != ""){
    query=paste(query,"AND sampleId in (",sampleString,")" )
  }
  
  raw_data = dbGetQuery(dbConnect, query)
}

get_tumorMarker_data<-function(dbConnect,sampleString="")
{
  query = paste("SELECT * FROM tumorMarker ",sep="")
  
  if (sampleString != ""){
    query=paste(query,"AND sampleId in (",sampleString,")" )
  }
  
  raw_data = dbGetQuery(dbConnect, query)
}

get_clinicalFindings_data<-function(dbConnect,sampleString="")
{
  query = paste("SELECT * FROM clinicalFindings ",sep="")
  
  if (sampleString != ""){
    query=paste(query,"AND sampleId in (",sampleString,")" )
  }
  
  raw_data = dbGetQuery(dbConnect, query)
}

get_treatmentResponse_data<-function(dbConnect,sampleString="")
{
  query = paste("SELECT * FROM treatmentResponse ",sep="")
  
  if (sampleString != ""){
    query=paste(query,"AND sampleId in (",sampleString,")" )
  }
  
  raw_data = dbGetQuery(dbConnect, query)
}

get_clinical_data<-function(dbConnect,sampleString="")
{
  query = paste("SELECT * FROM clinical ",sep="")
  
  if (sampleString != ""){
    query=paste(query,"AND sampleId in (",sampleString,")" )
  }
  
  raw_data = dbGetQuery(dbConnect, query)
}

get_treatment_data<-function(dbConnect,sampleString="")
{
  query = paste("select * from preTreatmentDrug inner join patient on patient.id = preTreatmentDrug.patientId ",sep="")
  
  if (sampleString != ""){
    query=paste(query,"AND sampleId in (",sampleString,")" )
  }
  
  raw_data = dbGetQuery(dbConnect, query)
}

get_gene_copyNumber<-function(dbConnect,sampleString="",geneString)
{
  query = paste("select sampleId,minCopyNumber from geneCopyNumber where minCopyNumber > -100",sep="")

  if (geneString != ""){
    query=paste(query," AND gene in ('",geneString,"')",sep="" )
  }  
  if (sampleString != ""){
    query=paste(query,"AND sampleId in (",sampleString,")" )
  }
  print(query)
  raw_data = dbGetQuery(dbConnect, query)
}



# Retrieve Data
load('~/hmf/RData/Processed/hpcDriversByGene.RData')
load('~/hmf/RData/Processed/highestPurityCohortSummary.RData')
dbConnect = dbConnect(MySQL(), dbname='hmfpatients_20180418', groups="RAnalysis")
clinicalData<-get_clinical_data(dbConnect)
treatmentData<-get_treatment_data(dbConnect)

clinicalFindingsData<-get_clinicalFindings_data(dbConnect)
# 3 types of clinical response
treatmentResponseData<-get_treatmentResponse_data(dbConnect)
tumorMarker<-get_tumorMarker_data(dbConnect)  # Bone metastasis measurements mainly for prostate cancers - 15% of patients. SHOULD Adjust to require a measurement exists
ranoData<-get_rano_data(dbConnect)  # Brain cancer measuments (eg. see RANO criteria for gliblastoma) - only in prod db - 1.5% of patients

HER2CopyNumber<-get_gene_copyNumber(dbConnect,,'ERBB2')

  

#immunotherapy  
# need to look more thoroughly - not just at first response but all treatment responses
# perhaps best response, perhaps time to first progession, etc.
# need to think carefully about biases with missing data
# neeed to check death date?
View(clinicalData %>% filter(treatmentType=='Immunotherapy') %>% group_by(treatment,firstResponse) %>% count() %>% spread(firstResponse,n,fill=0))

#Her2
# has the patient been treated with HER2 inhibitor?   Can we group combo chemo therapies?
View(clinicalData %>% filter(primaryTumorLocation=='Breast') %>% group_by(cancerSubtype,treatment) %>% 
       count() %>% spread(cancerSubtype,n,fill=0))
View(clinicalData %>% filter(primaryTumorLocation=='Breast',cancerSubtype %like% 'HER2-positive') %>% group_by(treatment,firstResponse) %>% count() %>% spread(firstResponse,n,fill=0))

View(hpcDriversByGene %>% filter(gene=='ERBB2',cancerType=='Breast'))
HER2positive = merge(clinicalData %>% filter(primaryTumorLocation=='Breast',cancerSubtype %like% 'HER2-positive'),
      hpcDriversByGene %>% filter(gene=='ERBB2',cancerType=='Breast'),
      by='sampleId',
      all.x=T)
View(merge(HER2positive,HER2CopyNumber,      
           by='sampleId',
           all.x=T))
View(HER2positive %>% group_by(firstResponse,driver) %>% count() %>% spread(firstResponse,n,fill=0))
View(HER2CopyNumber)


# Targeted - try loooking at VHL biomarker vs Pazopanib
View(clinicalData %>% filter(treatmentType=='Targeted therapy') %>% group_by(treatment,firstResponse) %>% count() %>% spread(firstResponse,n,fill=0))

