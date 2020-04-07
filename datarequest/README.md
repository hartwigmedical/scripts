
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
* [Other tips on working with GCP](#other-tips-on-working-with-GCP)

### General Notes

 - The data will be made available for GCP accounts of data download contacts. For these accounts multi-factor authentication needs to be enabled.
 - Please be aware that the sample IDs provided to you are not publicly shareable. In case you want to publish something with respect to a specific sample from the HMF database you need to use the HMF Sample ID for this. For newer DRs these have been made available; for older requests these can be requested for the samples you have been provided. 
 - Internally at HMF we load up all data into a MySQL database. The scheme and code to set this up yourself can be found on our [resources page](https://resources.hartwigmedicalfoundation.nl/).
 
Please use the **unique ID** given to your request (eg "DR-XXX") in any communication with us about your data request.

#### More information
- For source code of our analysis pipeline see our [pipeline5 repo](https://github.com/hartwigmedical/pipeline5).
- For source code of all HMF tools see our [hmftools repo](https://github.com/hartwigmedical/hmftools).
- For an example patient report see our [resources page](https://resources.hartwigmedicalfoundation.nl/).
- For various resource files used in the analysis see our [resources page](https://resources.hartwigmedicalfoundation.nl/).


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


### Accessing HMF Data

The process around accessing data is still evolving, but we now have a process to make all our BAMs available which was previously 
impossible. Given the size of the data, our key challenge is avoid any copying or duplication of the data. We do this by adding users 
directly to the ACL of each file they have access, and providing a manifest containing the URLs and parseable metadata they can use to 
download to a VM and organize appropriately.

The manifest is also still evolving, but in its current form it looks like:

```json
{
  "id": "example",
  "tertiary": [
    {
      "gsutilUrl": "gs://hmf-dr-example/clinical.tar.gz",
      "tags": {
        "moleculetype": "DNA",
        "datatype": "CLINICAL"
      }
    }
  ],
  "accounts": [
    {
      "email": "user@google.com"
    }
  ],
  "patients": [
    {
      "patientId": "COLO829v003",
      "samples": [
        {
          "sampleId": "COLO829v003T",
          "data": [
            {
              "gsutilUrl": "gs://hmf-output-test/COLO_VALIDATION_5_7_1314/COLO829v003R/aligner/COLO829v003R.bam",
              "tags": {
                "moleculetype": "DNA",
                "filetype": "BAM"
              }
            },
            {
              "gsutilUrl": "gs://hmf-output-test/COLO_VALIDATION_5_7_1314/COLO829v003R/aligner/COLO829v003R.bam.bai",
              "tags": {
                "moleculetype": "DNA",
                "filetype": "BAI"
              }
            },
            {
              "gsutilUrl": "gs://hmf-output-test/COLO_VALIDATION_5_7_1314/COLO829v003T/aligner/COLO829v003T.bam",
              "tags": {
                "moleculetype": "DNA",
                "filetype": "BAM"
              }
            },
            {
              "gsutilUrl": "gs://hmf-output-test/COLO_VALIDATION_5_7_1314/COLO829v003T/aligner/COLO829v003T.bam.bai",
              "tags": {
                "moleculetype": "DNA",
                "filetype": "BAI"
              }
            }
          ]
        }
      ]
    }
  ]
}
```

Downloading data outside of the `europe-west4` region carries a data egress cost of €0.11 per GB. We strongly recommend creating VMs close
to the data to avoid this cost and time of download. That said, it is possible to download, but you must specify your own project in the
`gsutil` command to do so. This is also necessary to download to a vm or another bucket, but you won't be charged anything if its in the
same region. For example:

```bash
gsutil -u my-project gs://hmf-output-test/COLO_VALIDATION_5_7_1314/COLO829v003R/aligner/COLO829v003R.bam ./
```

### Scaling Analysis

Running commands interactively can work for small workloads, but to get a real analysis done and really take advantage of GCP you'll want
to automate provisioning many VMs to scale up. The general pattern for this is:
* Create a VM with a predefined startup script.
* Within the startup script, download the data you need
* Within the startup script, run your analysis
* Within the startup scrupt, upload the results to your own bucket
* Terminate the VM

For this purpose, we're creating a tool you can use to take advantage of all our GCP automation, cost savings and perfomance improvements
called [Batch 5](https://github.com/hartwigmedical/pipeline5/blob/master/batch/README.md). Batch5 takes care of handling pre-emptible vms,
setting up local ssds, managing failures, exposing logs, etc.

We're also going to evaluate the Broad Institutes solution to this called [Terra](https://app.terra.bio/). Terra allows you to run pipelines
via common worflow languages (WDL, CWL), but also create notebooks. It runs on GCP and is a joint project between Broad and Google.



### Other tips working with GCP


#### Using the console

Your first entry point to GCP will likely be the console. It gives you a tool to both interact and monitor various GCP services and even
create a quick shell for more advances operations.

Hit the following link to start up the console: https://console.cloud.google.com/

### Projects, Users and Roles

When you begin in the console you'll notice that you are working within a *Project*. Project group together services and users along with 
billing. Most of the time you'll work in a single project along with close colleagues, and you're interactions with GCP will be confined
there.

Good to note that a user can belong to many projects.

IAM is used to organize users permissions. It uses the concepts of roles, which group permissions into cohesive abilities to perform tasks.

Unfortunately there are a lot of roles, and its not always clear how they map to a give permission. 
[This page](https://cloud.google.com/iam/docs/understanding-roles) gives the overview of how roles map to permissions.

### Using the command line

The fastest way to get to a shell is to use the "cloud shell" functionality. The cloud shell already has all tools installed and will 
default to the project you have open. It works by provisioning a small VM.

![Cloud Shell](https://github.com/hartwigmedical/gcpworkshop/blob/master/images/cloud-shell.png)

Open a cloud shell and try the following:

```
gcloud auth list
gsutil ls
```

That said, the shell still runs in the browser so can be a little slow and not suited to long running interactions. You can also do 
everything from a terminal on your own maching with the `gcloud` and `gsutil` tools.

Open a terminal and try the following:
```
gcloud auth login
gsutil ls
```

### Google Cloud Storage (GCS)

Storage in the cloud is one of its most power tools and concepts. GCS is an object storage system which works with buckets, paths and files.

The first time, its best to create a bucket via the console:

![Create a bucket](https://github.com/hartwigmedical/gcpworkshop/blob/master/images/storage-demo-1.png)

Important concepts:
* The region is where your bucket will be created. Best to pick a local region.
* The storage class impacts how your data is stored. Standard will be used for frequently accessed data. Nearline and Coldine for less 
accessed.
* The access control can be per bucket or per file, normally per bucket is fine.
* By default data in encrypted with Google's internal keys. Customer Managed Encryption is also available, but does carry a large 
performance penalty.

Create a bucket called `{yourname}-gcpdemo`. Bucket names have some restrictions: they must be globally unique; they must not have 
underscores; they must be lower-case. This is because you can expose data via direct URL access and must conform to the rules of DNS.

From there its best to move to the command line. After you've create your bucket, you can use the `gsutil` tool to upload, download, and
copy data between buckets. Try the following:

```
touch myfile.txt
gsutil -m cp myfile.txt gs://{yourname}-gcpdemo/
gsutil ls gs://{yourname}-gcpdemo/
```

The `-m` flag above will do the download in parallel. Keep in mind that this will be heavy on the bandwidth, but have a huge impact on the
speed of your operation. You can also get things to go even faster by uploading data as a composite object (aka multipart). The following
command will upload any file over 150M as a composite object, which also allows it to be downloaded as such. **Note**: You will not have an
MD5 checksum if a file is uploaded as composite.

```bash
gsutil -m -o GSUtil:parallel_composite_upload_threshold=150M cp myfile.txt gs://{yourname}-gcpdemo/
```

You can also browse your data from within the console. Find Storage in your console and you can see your new bucket and file there.

You can share your data with other users using the Access Control List of the bucket. Adding users is easiest via the console, but can also 
be done with `gsutil`. 


### Google Compute Engine (GCE)

GCE gives you the ability to create VMs, make images, and deploy containers. Next we'll create our own VM and learn how to access it via 
SSH.

Start with reviewing all the options from the Create VM on the [console](https://console.cloud.google.com/compute/instancesAdd):

![Create a VM](https://github.com/hartwigmedical/gcpworkshop/blob/master/images/compute-demo-1.png)

Important Concepts:
* VMs also need to deployed to a region, but also a zone, which is part of the datacenter
* There are many pre-defined machine types, all pricing is linear to CPU and memory. You can define custom machines via the command line.
* You can define a disk when you define your VM, but this can also be done independently, and read-only disks can be shared.
* Inside the VM, service interactions are performed by the project service account. A service account can run headless, but otherwise has
permissions just like a normal user. 

You can also create and manage vms via the command line. 

```bash
gcloud instances create my-instance --zone=europe-west4-a --boot-disk-size 100
```

The [VM instance overview](https://console.cloud.google.com/compute/instances) in the console is a nice way to see the status of all VMs:

![VM Instances](https://github.com/hartwigmedical/gcpworkshop/blob/master/images/compute-demo-2.png)

The SSH button gives you a quick browser shell to the machine. Again, the browser has its limitations, so sometimes it's better to tunnel 
into your VM from your own terminal:

```bash
gcloud compute --project "your-project" ssh --zone "europe-west4-a" "your-vm"
```

Try SSH'ing into your new VM and running the same `gsutil` commands we ran in the previous sections.

### Images

Images are a handy way to save and share state of a VM. You can take an image on the command line by:

```bash
 gcloud images create my-instance-image-1 --family=my-instance-image --source-disk=my-instance-image --source-disk-zone=europe-west4"
```

### GCE Cost Savings

During our migration to GCP, we found two major ways to cut back on our compute costs.

#### Pre-emptible VMS

The single biggest win we had was moving to all pre-emptible VMs. A pre-emptible VM, is one that GCP can claim back at any time, and will
always be shut-down at 24 hours. This is how Google sells excess capacity, and claims it back when they need it. It is a great fit for 
analysis workloads running between 1-24 hours, and the savings are big (80%). 

You can mark a VM pre-emptible on VM creation *Create VM -> Management -> Availability Policy* and through the `glcoud` CLI.

![Make a pre-emptible VM](https://github.com/hartwigmedical/gcpworkshop/blob/master/images/compute-demo-3.png)

While creating them is easy, its worth having some automation around them to manage pre-emptions. In our pipline, we handle the pre-empted
signal by polling the GCE API, then re-starting the workload in a new zone.

#### Local SSDs

Smaller impact and more complex, local SSDs give you the best disk performance at the lowest cost. That said, you need to perform some steps
if you need a disk larger than 375GB and data will be completely transient. A VM with local SSDs cannot be restarted and all data there is
lost. 

They are a great way to both speed up and reduce cost of a completely transient workload.


