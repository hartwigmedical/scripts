
HMF Data Request Guide
===== 

This page provides information about the data you received within the context of a data request (DR) from Hartwig Medical Foundation.

### Contents

* [General notes](#general-notes)
* [Sample selection](#sample-selection)
* [Tips on downloading files ](#tips-on-downloading-files)

### General Notes
 - Sharing of data is done via an OKTA account which requires dual-factor authentication. The OKTA account is used to access our [Nextcloud Portal](https://nc.hartwigmedicalfoundation.nl) and our [Download Portal](https://portal.hartwigmedicalfoundation.nl).
 - Please be aware that the sample IDs provided to you are not publicly shareable. In case you want to publish something with respect to a specific sample from the HMF database you need to use the HMF Sample ID for this. These can be requested for the samples you have been provided. 
 - Internally at HMF we load up all data into a MySQL database. The scheme and code to set this up yourself can be found on our [resources page](http://resources.hartwigmedicalfoundation.nl).
 - Methods used to collect and analyse our samples can be found on our [methods](./README_METHODS.md) page.
 
Please use the **unique ID** given to your request (eg "DR-XXX") in any communication with us about your data request.

### Clinical Data

Clinical data is shared in a **metadata.tar** via [Nextcloud Portal](https://nc.hartwigmedicalfoundation.nl).

Some notes about the clinical data:
- As much information as possible from the Electronic Case Report Form (eCRF) of the respective clinical studies is gathered. Please be aware that records are not guaranteed to be complete.
- For patients who participated in the DRUP study we can not share any treatment related information.
- The "purity" field in the metadata is the percentage tumor cells derived from WGS data by [PURPLE](https://github.com/hartwigmedical/hmftools/tree/master/purity-ploidy-estimator).

### Somatic Data

Somatic data is shared via a url to **somatics.tar** via [Nextcloud Portal](https://nc.hartwigmedicalfoundation.nl).

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

### Germline Data

Germline data is shared via a url to **germline.tar** via [Nextcloud Portal](https://nc.hartwigmedicalfoundation.nl).

We share the SNVs and small INDELs called from the reference sample using GATK haplotype caller.

### Alignments (BAM files)

Sharing of BAM files is currently only supported for samples which have been previously published, but in case you do have access they can be accessed in the [Download Portal](https://portal.hartwigmedicalfoundation.nl).

**Some notes to keep in mind:**
- The session time is 2 hours, after this period you will have to login again
- Files can be very large (up to 300GB) so please check before starting a download!
- When you create links to files, be aware that these links are valid for 24 hours. After this period you will have to create new links.

**Example loading BAM file in IGV:**
- Create the links for the BAM and accompanying BAI file (by clicking the most right icon next to a run)
- Open the IGV program and choose option "Load from URL"
- Paste the BAM file link at field "File URL"
- Paste the BAI file link at field "Index URL"

### Sample selection

By default, in addition to data-request specific criteria, samples for which one of the below applies are **excluded**:

- Samples from patients where informed consent is from before 21 April 2016.
- Samples with poor quality (PURPLE qcStatus != PASS).
- Samples without any tumor evidence (PURPLE status = NO_TUMOR).
- Samples with less than 19.5% tumor cells (PURPLE purity < 0.195).

### Tips on downloading files 

Links to URLs can be retrieved or created via [Nextcloud Portal](https://nc.hartwigmedicalfoundation.nl) and [Download Portal](https://portal.hartwigmedicalfoundation.nl). For smaller files, wget generally works fine and is provided on standard linux-based operating systems.

**Example download using WGET:**

- Copy the links into a new text file (eg links.txt)
- Use the following command to download: 

```sh
wget --content-disposition -i links.txt
```

For bigger files or if wget doesn't work in the first place, we recommend to use the [aria2 download tool](https://aria2.github.io/). This will allow you to download with multiple streams speeding up the download significantly.

**Example download using aria2 (with 5 connections):**
```sh
aria2c -x 5 -c -i links.txt
```

### More information
- For source code of our analysis pipeline see our [pipeline5 repo](https://github.com/hartwigmedical/pipeline5).
- For source code of all HMF tools see our [hmftools repo](https://github.com/hartwigmedical/hmftools).
- For an example patient report see our [resources page](https://resources.hartwigmedicalfoundation.nl/).
- For various resource files used in the analysis see our [resources page](https://resources.hartwigmedicalfoundation.nl/).
