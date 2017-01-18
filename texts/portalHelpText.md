## Hartwig Medical Foundation Download Portal

**The portal**

Via this portal you have access to molecular data from your group or center. For each sample/patient certain files can be downloaded by creating links (valid for 24 hours). The links can be used to download (via browser or standard tools like WGET/cURL) or, in the case of indexed files like the BAM format, directly read ("stream") into a program like IGV (Integrative Genome Browser) so you do not have to download large files completely just to check a small part of them.

[https://portal.hartwigmedicalfoundation.nl/](https://portal.hartwigmedicalfoundation.nl/)

**Good to know**

- Files can be very large (up to 300GB) so please check before starting a download!
- The session time is 2 hours, after this period you will have to login again
- Links are valid for 24 hours, after this period you will have to create new links
- This portal contains "molecular" data, no clinical information (for data requests that require clinical information please see our [Data Policy](https://www.hartwigmedicalfoundation.nl/databeleid/))
- Examples of available data formats are BAM, gVCF and VCF
- See below for some usage examples

Please let us know if you have comments or want to share your portal experience. For more information and contact please visit our [website](https://www.hartwigmedicalfoundation.nl).

-----
### Usage examples

**Example loading BAM file in IGV:**
- create the links for the BAM and accompanying BAI file (by clicking the most right icon next to a run)
- open the IGV program and choose option "Load from URL"
- paste the BAM file link at field "File URL"
- paste the BAI file link at field "Index URL"

Note: same procedure applies to VCF files

**Example download using WGET:**
- create the links (by clicking the most right icon next to a run)
- copy the links into a new text file (eg download.txt)
- use the following command to download: 

```sh
wget --content-disposition -i download.txt
```

Note: we also recommend to check out the [aria2 download tool](https://aria2.github.io/). This will allow you to download with mutiple streams speeding up the download significantly.

