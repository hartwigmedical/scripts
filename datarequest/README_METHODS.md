HMF Description of Methods
===== 

This page describes the methods used to generate the data serviced to researchers through data requests (DRs).

## Sample collection

Patients with advanced cancer not curable by local treatment options and being candidates for any type of systemic treatment and any line of treatment were included as part of the CPCT-02 (NCT01855477) and DRUP (NCT02925234) clinical studies, which were approved by the medical ethical committees (METC) of the University Medical Center Utrecht and the Netherlands Cancer Institute, respectively. 

A total of 41 academic, teaching and general hospitals across the Netherlands participated in these studies and collected material and clinical data by standardized protocols. Patients have given explicit consent for whole genome sequencing and data sharing for cancer research purposes. Clinical data, including primary tumor type, biopsy location, gender and birth year were collected in electronic case record forms and stored in a central database.

Core needle biopsies were sampled from the metastatic lesion, or when considered not feasible or not safe, from the primary tumor site when still in situ. One to four biopsies were collected (average of 2.1 per patient) and frozen in liquid nitrogen directly after sampling. In parallel, a tube of blood was collected in CellSave or Streck tubes, which was shipped by room temperature to the central sequencing facility at the Hartwig Medical Foundation. Left-over material (biopsy, DNA) after sample processing was stored in biobanks associated with the studies at the University Medical Center Utrecht and the Netherlands Cancer Institute.

## Sequencing workflow

DNA was isolated from biopsy and blood on an automated setup (QiaSymphony) according to supplier's protocols (Qiagen) using the DSP DNA Midi kit for blood and QIAsymphony DSP DNA Mini kit for tissue and quantified (Qubit). Before starting DNA isolation from tissue, the biopsy was dissolved in 100 microliter Nuclease-free water by using the Qiagen TissueLyzer and split in two equal fractions for parallel automated DNA and RNA isolation (QiaSymphony). 

Typically, DNA yield for the tissue biopsy ranged between 50 and 5,000 ng. A total of 50 - 200 ng of DNA was used as input for TruSeq Nano LT library preparation (Illumina), which was performed on an automated liquid handling platform (Beckman Coulter). DNA was sheared using sonication (Covaris) to average fragment lengths of 450 nt. Barcoded libraries were sequenced as pools on HiSeq X (V2.5 reagents) and Novaseq 6000 S4 Reagent Kit generating 2 x 151 read pairs using standard settings (Illumina).

## Bioinformatics workflow

 BCL output from the HiSeqX and Novaseq6000 platform was converted using bcl2fastq tool (Illumina, versions 2.17 to 2.20 have been used) using default parameters. 

Reads were mapped to the reference genome GRCh37 using BWA-mem v0.7.x. Duplicates were marked for filtering. In addition, for all HiSeq X data, INDELs were realigned and base qualities were recalibrated prior to somatic SNV/INDEL calling using GATK v3.x. These steps have also been applied to some Novaseq data.  

GATK HaplotypeCaller v3.x was run to call germline variants in the reference sample, after a which a set of soft filters are applied which can be found in the germline VCFs of the samples. 

Strelka v1 was run to call somatic SNVs and small INDELs using non-default (reduced) filters to increase sensitivity. In addition various additional filters are run post-strelka to increase precision.

GRIDSS v2.x was run to call structural variants with various additional annotations applied. 

PURPLE v2.x was run to fit purity and ploidy.    

For more detailed information on tools and parameters used:
 - For germline SNV/INDEL calling and somatic SNV/INDEL calling please see the 'Methods' section of [HMF Pan Cancer Paper](https://www.nature.com/articles/s41586-019-1689-y)
 - For GRIDSS/Purple/LINX see [HMF Toolkit Paper](https://www.biorxiv.org/content/10.1101/781013v1)