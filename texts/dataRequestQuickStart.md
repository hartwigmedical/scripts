
Data Request Quick Start
===== 

Here you find information about the data you received within a data request from Hartwig Medical Foundation. Please use the unique ID given to your request (eg "DR-XXX") in any communication.

<!--
Most recent version of this text can be found at one the following repositories:  
[https://github.com/hartwigmedical/scripts/tree/master/texts](https://github.com/hartwigmedical/scripts/tree/master/texts)  
[https://github.com/hartwigmedical/texts/](https://github.com/hartwigmedical/texts/)
-->

In case your data request involves somatic data you can find these in our nextcloud portal:  
[https://nc.hartwigmedicalfoundation.nl/](https://nc.hartwigmedicalfoundation.nl/)

In case your data request involves germline data, you can find the Germline VCFs in our download-portal:   
[https://portal.hartwigmedicalfoundation.nl/](https://portal.hartwigmedicalfoundation.nl/)

Note: both download-portal and nextcloud are behind the same dual factor login (that works with app "OKTA Verify").

### Somatic data

Somatic data is shared in a gzipped tar via our [Nextcloud](https://nc.hartwigmedicalfoundation.nl/).

##### Contents of the gzipped tar file:
- ./metadata directory (content depends on the specifics of the request).
- ./data directory with various somatic data files per biopsy/set.

##### Per biopsy the following files are present:
- post_processed.vcf.gz (somatic SNVs and small INDELs).
- purple.sv.ann.vcf.gz (somatic structural variants).
- purple.cnv (somatic copy number regions).
- purple.gene.cnv (somatic copy number per gene).
- purple.purity (implied percentage of tumor cells).
- purple.purity.range (in depth information about purity measure).
- purple.qc (quality control outcome).
- purple.version (purple version used).
- circos.png (genome wide plot).
- [in case germline level is part of request] purple.germline.cnv (germline copy number regions).


##### Notes about the (clinical) metadata:
- We gather as much information as possible from the Electronic Case Report Form (eCRF) of the respective clinical studies, but please be aware that records are by no means complete.
- For biopsies from the DRUP study the "patientId" is replaced with the respective CPCT-patientId in case the patient was already known within the CPCT02 study (this way it is clear that these are from the same individual).
- For patients from the DRUP study we can not share any treatment related information.
- The "purity" mentioned is the percentage tumor cells derived from WGS data by the tool [PURPLE](https://github.com/hartwigmedical/hmftools/tree/master/purity-ploidy-estimator).


### Germline data

In case Germline level is included in your request, Germline variant VCFs (GATK based) are available via our [download-portal](https://portal.hartwigmedicalfoundation.nl/).

- On top of the download-portal interface you can find some download tips behind the question mark icon. Best is to use the aria2 config option as it allows you to select multiple sets, download a plain text config file and then use the aria2 download tool to perform download and md5sum file integrity check in one go.
- In case you want to use other command line tools like curl/wget to download files, make sure to put quotes around the URLs as certain characters might otherwise break the URL.

Sharing of BAM files is currently only supported for small numbers, but in case you do have access they are also visible in the download-portal.


### Data inclusion/exclusion

In addition to data-request specific criteria, by default we exclude datasets for which one of below applies.

- We exclude biopsies from patients where informed consent is from before 21 April 2016.
- We exclude biopsies with poor quality (PURPLE qcStatus != PASS).
- We exclude biopsies without any tumor evidence (PURPLE status = NO_TUMOR).
- We exclude biopsies with less than 20% tumor cells (PURPLE purity < 0.2).

### More detailed information
- Source code of our analysis pipeline: [https://github.com/hartwigmedical/pipeline](https://github.com/hartwigmedical/pipeline)
- Source code of all HMF tools: [https://github.com/hartwigmedical/hmftools](https://github.com/hartwigmedical/hmftools)
- For an explanation of most output see [PURPLE](https://github.com/hartwigmedical/hmftools/tree/master/purity-ploidy-estimator)
- For an example patient report see our [resources page](http://resources.hartwigmedicalfoundation.nl/)
- For various resource files used in our analysis see our [resources page](http://resources.hartwigmedicalfoundation.nl/)

### Final Notes
- We internally use a MySQL database and can share the scheme with you to set this up yourself as well. Please contact one of us for instructions. There is an setup example available at our [resources page](http://resources.hartwigmedicalfoundation.nl/).
