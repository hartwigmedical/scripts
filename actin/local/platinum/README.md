## Running COLO pipeline via platinum

Steps:
1. Clone the platinum repo from https://github.com/hartwigmedical/platinum, and build locally (see README)
2. Connect to VPN
3. Edit the [COLO829.yaml](COLO829.yaml) file with respect to the following fields:
   - image: Take the pipeline5 release you want to use (based on a tag)
   - outputBucket: Pick a unique, non-existing bucket where output will be written to
   - imageName: This refers to the VM image that will be used by your pipeline image.
4. From platinum repo, run "./platinum run -i /path/to/COLO829.yaml -n ${run_id}"
   - Note: ${run_id} should be short, eg "kd-240307-colo"

You can use the same method to run any other sample (either from BAM or from FASTQ) with any pipeline version, assuming you know where
the input data is stored, and `hmf-crunch` has access to it. 