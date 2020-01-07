HMF Description of Methods
===== 

This page describes the methods used to generate the data serviced to researchers through data requests (DRs). The samples in the Hartwig database are collected via three clinical studies: CPCT-02, DRUP-01 and WIDE.
Patients in these studies have given explicit consent for whole genome sequencing and data sharing for cancer research purposes.

These studies deliver clinical data to Hartwig to accompany the whole genome analysis. The clinical data is curated and standardized by Hartwig as part of the data ingestion algorithms found in our [Patient DB](https://github.com/hartwigmedical/hmftools/tree/master/patient-db)   

## CPCT-02 ([NCT01855477](https://clinicaltrials.gov/ct2/show/NCT01855477))

Patients with advanced cancer not curable by local treatment options and being candidates for any type of systemic treatment and any line of treatment were included as part of the Center for Personalized Cancer Treatment (CPCT-02) clinical study, which was approved by the medical ethical committee (METC) of the University Medical Center Utrecht. 

A total of 44 academic, teaching and general hospitals across the Netherlands participated in this study and collected material and clinical data by standardized protocols. 

Core needle biopsies were sampled from the metastatic lesion, or when considered not feasible or not safe, from the primary tumor site when still in situ. One to four biopsies were collected (average of 2.1 per patient) and frozen in liquid nitrogen directly after sampling. In parallel, a tube of blood was collected in CellSave or Streck tubes, which was shipped by room temperature to the central sequencing facility at the Hartwig Medical Foundation. Left-over material (biopsy, DNA) after sample processing was stored in a biobank at the University Medical Center Utrecht.

#### CPCT-02 clinical data collection

CPCT clinical data is collected in an electronic Case Record Form by local datamanagers at the 44 different hospitals, and supervised by central datamanagement. The current data extract has not been formally cleaned by CPCT. They are in the process of checking for any inconsistencies in the data and working with the hospitals to resolve these issues. Therefore, there may be inaccuracies in the data that are not resolved yet. Please carefully consider this when interpreting the data.

CPCT's standard quality control system consists of source data verification. This means about 10% of all data will be verified against source data at the different hospitals. There is an average of 10-20% of items that are not consistent with source documentation. This means the data will not be fully complete. In case of insufficient documentation by the treating physician there is no way to correct for this. 

Central datamanagement is in the process of performing a logic check, meaning that after data entry the consistency of the data is checked (eg are the dates possible and logical, is the data entered conform protocol) and if possible errors are resolved. The last check will be a medical check to ensure the data from each patient makes sense from the medical point of view. CPCT expects to have these checks completed for the entire dataset in the beginning of 2020.

There are specific considerations to keep in mind with different types of data:

 - Tumor type: The data on tumor type is reliable and will not contain more than 1-2% errors. CPCT has categorized tumors by tissue type of origin and provides one layer of subtype classification. Further histological subclassification is typically not present as the biopsies were collected from metastases and not part of routine diagnostics.
 
 - Prior treatment: The type of prior treatment will in general be reliable and whether the prior treatment fits with the diagnosis will be part of the medical check of the data. Best response to treatment is a retrospective feature and lacks full consistency. CPCT does not have the means to fully check this data. The dates that are given in many cases are estimated dates. Please be careful in interpreting the data on response and exposure times. If you want to make correlations between genomic data and duration of exposure / response these data are not completely verifiable.

 - Post biopsy treatment: From here the data is collected prospectively and thus more reliable and uniformly present.

 - Response to treatment after biopsy: These data are collected prospectively and should be accurate. There is however the need to understand that CPCT did not mandate the exact time at which the evaluation was performed. Therefore differences in parameters such as time to progression are not a reliable measure of efficacy (as the physician choice to perform an evaluation biases the data). Also CPCT did not mandate confirmation of responses so there will be an underestimation of the number of confirmed responses. The time on treatment is a reliable measure, however incorporates the stopping of medication due to toxicity. CPCT does not make a distinction between these 2 stop reasons. Also note that CPCT did not collect the data on dose intensity.

 - Survival: the data on survival are not present for every patient. Connecting the data to survival data is something CPCT aims to be performing in 2020.

## DRUP-01 ([NCT02925234](https://clinicaltrials.gov/ct2/show/NCT02925234))

Patients who were eligible for the study had an advanced or metastatic solid tumour, multiple myeloma or B cell non-Hodgkin lymphoma, and had exhausted standard-treatment options. A tumour genetic or protein-expression test (CPCT or regular diagnostics) had to have revealed a potentially actionable variant, for which FDA- and/or EMA-approved targeted (or immuno-) therapy was available—but not for the tumour type in question. This clinical study was approved by the Medical Ethical Committee of the Netherlands Cancer Institute in Amsterdam

Upon case submission, the study team attempted to match each patient to the appropriate study treatment, according to pre-defined matching rules. If a slot for a matching study treatment was available (to which the patient consented) the patient could be enrolled, if all drug-specific selection criteria were met. Once a variant-drug match was identified, all patients underwent a fresh baseline tumor biopsy for WGS. 

Core needle biopsies were sampled similar according to the CPCT protocol as outlined above and also shipped to the Hartwig Medical Foundation for sequencing. Left-over material (biopsy, DNA) after sample processing was stored in a biobank at the Netherlands Cancer Institute in Amsterdam. 

A total of 37 academic, teaching and general hospitals across the Netherlands participated in this study and collected material and clinical data by standardized protocols.
DRUP clinical data collection

#### DRUP clinical data collection

DRUP clinical data is collected in an electronic Case Record Form by local datamanagers at the 37 different hospitals, and supervised by central datamanagement. The current data extract has not been fully cleaned by DRUP. They are in the process of checking for any inconsistencies in the data and working with the hospitals to resolve these issues. Therefore, there may be inaccuracies in the data that are not resolved yet. Please carefully consider this when interpreting the data.

DRUP’s standard quality control system consists of source data verification. This means that every first patient for a certain drug in a hospital is verified again source data and then about 10% randomly after. In case of insufficient documentation by the treating physician there is no way to correct for this.

Central datamanagement is in the process of performing a logic check, meaning that after data entry the consistency of the data is checked (e.g. are the dates possible and logical, is the data entered conform protocol) and if possible errors are resolved. The last check will be a medical check to ensure the data from each patient makes sense from the medical point of view. 

There are specific considerations to keep in mind with different types of data:

- Tumor type: The data on tumor type is reliable and will not contain more than 1-2% errors. DRUP has categorized tumors by tissue type of origin and provides one layer of subtype classification. Further histological sub classification is typically not present as the biopsies were generally collected from metastases and not part of routine diagnostics.

- Prior treatment: The type of prior treatment will in general be reliable and whether the prior treatment fits with the diagnosis will be part of the medical check of the data. Best response to treatment is a retrospective feature and lacks full consistency. DRUP does not have the means to fully check this data. The dates that are given in many cases are estimated dates. Please be careful in interpreting the data on response and exposure times. If you want to make correlations between genomic data and duration of exposure / response these data are not completely verifiable.

- Treatment data: the DRUP trial is an ongoing study. Therefore, clinical, treatment outcome and safety data cannot (yet) be issued.  

- Survival: the data on survival are not present for every patient. 


## WIDE ()

Coming soon...

## Sequencing workflow

DNA was isolated from biopsy and blood on an automated setup (QiaSymphony) according to supplier's protocols (Qiagen) using the DSP DNA Midi kit for blood and QIAsymphony DSP DNA Mini kit for tissue and quantified (Qubit). Before starting DNA isolation from tissue, the biopsy was dissolved in 100 microliter Nuclease-free water by using the Qiagen TissueLyzer and split in two equal fractions for parallel automated DNA and RNA isolation (QiaSymphony). 

Typically, DNA yield for the tissue biopsy ranged between 50 and 5,000 ng. A total of 50 - 200 ng of DNA was used as input for TruSeq Nano LT library preparation (Illumina), which was performed on an automated liquid handling platform (Beckman Coulter). DNA was sheared using sonication (Covaris) to average fragment lengths of 450 nt. Barcoded libraries were sequenced as pools on HiSeq X (V2.5 reagents) and Novaseq 6000 S4 Reagent Kit generating 2 x 151 read pairs using standard settings (Illumina).

## Bioinformatics workflow

 BCL output from the HiSeqX and Novaseq6000 platform was converted using bcl2fastq tool (Illumina, versions 2.17 to 2.20 have been used) using default parameters. 

Reads were mapped to the reference genome GRCh37 using BWA-mem v0.7.x. Duplicates were marked for filtering. In addition, for all HiSeq X data, INDELs were realigned and base qualities were recalibrated prior to somatic SNV/INDEL calling using GATK v3.x. These steps have also been applied to some Novaseq data.  

GATK HaplotypeCaller v3.x was run to call germline variants in the reference sample, after a which a set of soft filters are applied which can be found in the germline VCFs of the samples. 

Strelka v1 was run to call somatic SNVs and small INDELs using non-default (reduced) filters to increase sensitivity. In addition various additional filters are run post-strelka to increase precision.

[GRIDSS](https://github.com/PapenfussLab/gridss) v2.x was run to call structural variants with various additional annotations applied. 
  
[PURPLE](https://github.com/hartwigmedical/hmftools/tree/master/purity-ploidy-estimator) v2.x was run to fit purity and ploidy.
  
[LINX](https://github.com/hartwigmedical/hmftools/tree/master/sv-linx) v1.x was run to analyse the structural variants.    
  
For more detailed information on tools and parameters used:
 - For germline SNV/INDEL calling and somatic SNV/INDEL calling please see the 'Methods' section of [HMF Pan Cancer Paper](https://www.nature.com/articles/s41586-019-1689-y)
 - For GRIDSS/Purple/LINX see [HMF Toolkit Paper](https://www.biorxiv.org/content/10.1101/781013v1)