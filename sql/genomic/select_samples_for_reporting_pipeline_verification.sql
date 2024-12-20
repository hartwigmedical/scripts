# HRD
select * from chord inner join purity on purity.sampleId = chord.sampleId where hrStatus = "CANNOT_BE_DETERMINED" order by modified desc;
select * from chord inner join purity on purity.sampleId = chord.sampleId where hrStatus = "HR_DEFICIENT" order by modified desc;
select * from chord inner join purity on purity.sampleId = chord.sampleId where hrStatus = "HR_PROFICIENT" order by modified desc;

#virus
select * from virusAnnotation where reported and interpretation = "HPV" and likelihood = "HIGH" order by modified desc;

# MSI
select * from purity where msStatus = "MSS" order by modified desc;
select * from purity where msStatus = "MSI" order by modified desc;

#TMB
select * from purity where tmbStatus = "LOW" and qcStatus = "WARN_LOW_PURITY" order by modified desc;
select * from purity where tmbStatus = "HIGH" and qcStatus = "PASS" order by modified desc;

# fusion
select * from svFusion where reported = 1 and reportedType = "KNOWN_PAIR" order by modified desc;
select * from svFusion where reported = 1 and reportedType like '%PROMISCUOUS%' order by modified desc;
select * from svFusion where reported = 1 and reportedType = "EXON_DEL_DUP" order by modified desc;

#mutation
select * from driverCatalog where driver = "MUTATION" order by modified desc;
select * from driverCatalog where driver = "GERMLINE_MUTATION" order by modified desc;

#SV
select * from driverCatalog where driver = "DEL" order by modified desc;
select * from driverCatalog where driver = "AMP" order by modified desc;
select * from driverCatalog where driver = "PARTIAL_AMP" order by modified desc;
select * from driverCatalog where driver = "DISRUPTION" order by modified desc;
select * from driverCatalog where driver = "HOM_DEL_DISRUPTION" order by modified desc;
select * from driverCatalog where driver = "GERMLINE_DISRUPTION" order by modified desc;
select * from driverCatalog where driver = "GERMLINE_DELETION" order by modified desc;
select * from driverCatalog where driver = "HOM_DUP_DISRUPTION" order by modified desc;