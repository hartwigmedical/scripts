
Quick Start Data Request
===== 

This page provides information about the data you received within a data request from Hartwig Medical Foundation. Please use the **unique ID** given to your request (eg "DR-XXX") in any communication.

<!--
Most recent version of this text can be found at one the following repositories:  
[https://github.com/hartwigmedical/scripts/tree/master/texts](https://github.com/hartwigmedical/scripts/tree/master/texts)  
[https://github.com/hartwigmedical/texts/](https://github.com/hartwigmedical/texts/)
-->

In case your data request involves tertiary analysis data (VCF/TXT) you can find these in our [Nextcloud portal](https://nc.hartwigmedicalfoundation.nl/). In case your data request involves BAM files, you can find the links to the BAMs in our [Download-portal](https://portal.hartwigmedicalfoundation.nl/). Both the download-portal and nextcloud are behind the same dual factor login (with app "OKTA Verify").

### Somatic

Somatic data is shared in a **gzipped tar** via our [Nextcloud](https://nc.hartwigmedicalfoundation.nl/).

##### Contents of the tar file:
- ./metadata directory (content depends on the specifics of the request).
- ./somatics directory with various somatic data files per biopsy.
- ./germline directory with germline calls of the reference blood sample. *[in case germline level is part of request]*

##### Per biopsy the following files are present:
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

##### Notes about the (clinical) metadata:
- As much information as possible from the Electronic Case Report Form (eCRF) of the respective clinical studies is gathered. Please be aware that records are not guaranteed to be complete.
- For patients from the DRUP study we can not share any treatment related information.
- The mentioned "purity" is the percentage tumor cells derived from WGS data by [PURPLE](https://github.com/hartwigmedical/hmftools/tree/master/purity-ploidy-estimator).


### Germline

In case Germline level is included in your request, GATK based variant calls of the healthy reference sample (blood) are available via our [Nextcloud](https://nc.hartwigmedicalfoundation.nl/).

### Alignments

Sharing of BAM files is currently only supported for small numbers, but in case you do have access they are also visible in the download-portal.

### Sample selection

By default, in addition to data-request specific criteria, samples for which one of the below applies are **exluded**:

- Samples from patients where informed consent is from before 21 April 2016.
- Samples with poor quality (PURPLE qcStatus != PASS).
- Samples without any tumor evidence (PURPLE status = NO_TUMOR).
- Samples with less than 19.5% tumor cells (PURPLE purity < 0.195).

### More information
- For source code of our analysis pipeline see our [pipeline5 repo](https://github.com/hartwigmedical/pipeline5).
- For source code of all HMF tools see our [hmftools repo](https://github.com/hartwigmedical/hmftools).
- For an explanation of most output see our [PURPLE tool](https://github.com/hartwigmedical/hmftools/tree/master/purity-ploidy-estimator).
- For an example patient report see our [resources page](http://resources.hartwigmedicalfoundation.nl/).
- For various resource files used in the analysis see our [resources page](http://resources.hartwigmedicalfoundation.nl/).

### Final Notes
- At Hartwig Medical Foundation we use a MySQL database, the scheme and code to set this up yourself can be shared. Please contact one of us for instructions. There is small example available at our [resources page](http://resources.hartwigmedicalfoundation.nl/).
