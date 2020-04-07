
HMF Data Request Guide
===== 

This page provides practical information on how to access / work with the data that will be made available to you within the context of a Data Access Request (DR) from Hartwig Medical Foundation. 

Hartwig Medical foundation will make the licenced data available through the Google Cloud Platform (GCP). 

Note: more details on the methods used to generate both the genomic and clinical data can be found on a separate [Methods](./README_METHODS.md) page.


### Contents

* [General notes](#general-notes)
* [Sample selection](#sample-selection)
* [Format of the data made available](#format-of-the-data-made-available)
* [Tips on accessing the data trough GCP](#tips-on-accessing-the-data-trough-GCP)


### General Notes

 - The data will be made available for GCP accounts of data download contacts. For these accounts multi-factor authentication needs to be enabled.
 - Please be aware that the sample IDs provided to you are not publicly shareable. In case you want to publish something with respect to a specific sample from the HMF database you need to use the HMF Sample ID for this. For newer DRs these have been made available; for older requests these can be requested for the samples you have been provided. 
 - Internally at HMF we load up all data into a MySQL database. The scheme and code to set this up yourself can be found on our [resources page](https://resources.hartwigmedicalfoundation.nl/).
 
Please use the **unique ID** given to your request (eg "DR-XXX") in any communication with us about your data request.


### Sample selection

By default, in addition to data-request specific criteria, samples for which one of the below applies are **excluded**:

- Samples from patients where informed consent is from before 21 April 2016.
- Samples with poor quality (PURPLE qcStatus != PASS).
- Samples without any tumor evidence (PURPLE status = NO_TUMOR).
- Samples with less than 19.5% tumor cells (PURPLE purity < 0.195).



### Format of the data made available

#### Clinical Data (TSV format)

Clinical data will be made available in a **metadata.tar** via GCP.

Some notes about the clinical data:
- As much information as possible from the Electronic Case Report Form (eCRF) of the respective clinical studies is gathered. Please be aware that records are not guaranteed to be complete.
- For patients who participated in the DRUP study we can not share any treatment related information.
- The "purity" field in the metadata is the percentage tumor cells derived from WGS data by [PURPLE](https://github.com/hartwigmedical/hmftools/tree/master/purity-ploidy-estimator).

Please find more details on the methods used to generate both the genomic and clinical data on a separate [Methods](./README_METHODS.md) page.

#### Somatic Data (VCF/TXT formats)

Somatic data will be made available in a **somatics.tar** via GCP.

##### Per sample the following files are present:
- purple.somatic.vcf.gz (somatic SNVs and small INDELs).
- purple.sv.ann.vcf.gz (somatic structural variants).
- purple.cnv.somatic.tsv (somatic copy number regions).
- purple.cnv.gene.tsv (somatic copy number per gene).
- purple.purity.tsv (implied percentage of tumor cells).
- purple.purity.range.tsv (in depth information about purity measure).
- purple.qc (quality control outcome).
- purple.version (purple version used).
- circos.png (genome wide plot).
- driver.catalog.tsv (affected driver genes).
- purple.cnv.germline.tsv (germline copy number regions). *[in case germline level is part of request]*

For an explanation of the contents of these files, see [PURPLE](https://github.com/hartwigmedical/hmftools/tree/master/purity-ploidy-estimator).  

#### Germline Data (VCF/TXT formats)

Germline data will be made available in a **germline.tar** via GCP.

We share the SNVs and small INDELs called from the reference sample using GATK haplotype caller.

### Aligned readout data (CRAM format)

Aligned readout data will be made available per sample via GCP.

**Some notes to keep in mind:**
- Files can be very large (up to 300GB) so consider this when you want to egress the files (comes with costs).

**Example loading BAM file in IGV:**
- Create the links for the CRAM and accompanying CRAI file (by clicking the most right icon next to a run)
- Open the IGV program and choose option "Load from URL"
- Paste the BAM file link at field "File URL"
- Paste the BAI file link at field "Index URL"

Please find more details on the methods used to generate both the genomic and clinical data on a separate [Methods](./README_METHODS.md) page.

#### RNAseq data (FASTQ format)

RNAseq data will be made available per sample via GCP.


### Tips on accessing the data trough GCP


**Example egress data from GCP:**
```sh

```

### More information
- For source code of our analysis pipeline see our [pipeline5 repo](https://github.com/hartwigmedical/pipeline5).
- For source code of all HMF tools see our [hmftools repo](https://github.com/hartwigmedical/hmftools).
- For an example patient report see our [resources page](https://resources.hartwigmedicalfoundation.nl/).
- For various resource files used in the analysis see our [resources page](https://resources.hartwigmedicalfoundation.nl/).
