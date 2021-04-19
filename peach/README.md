# PEACH

**P**harmacogenetic **E**valuator **A**nd **C**aller of **H**aplotypes (PEACH) is a pharmacogenomics tool developed for the [Hartwig Medical Foundation pipeline](https://github.com/hartwigmedical/pipeline5).
It imports haplotypes and related variants from a curated JSON file, and inspects the presence of these variants in a germline VCF. 

It creates two output files:
* [sample].peach.genotype.tsv; contains on each line a determined genotype of the sample for a specific gene, expressed in terms of haplotypes.
* [sample].peach.calls.tsv; contains all the variants from the JSON file and their respective calls and filters.
 
## Contents

* [Installation](#installation)
* [Arguments](#arguments)
  + [Mandatory Arguments](#mandatory-arguments)
  + [Optional Arguments](#optional-arguments)

## Installation
If you want to run PEACH, please generate a local Python 3 venv and install the requirements:

```bash
$ python3 -m venv [path/to/new/virtual/environment, for example: ./peach]
$ source [path/to/new/venv, for example: ./peach/bin/activate]
(peach) $ pip install -r requirements.txt
```

## Arguments
Remember to source the virtualenv before running `main.py`.

####Example Usage
```
(peach) $ python main.py \
    input.vcf.gz \
    COLO829T \
    COLO829R \
    1.0 \
    /path/to/outputdir/ \
    /path/to/panel.json \
    /path/to/vcftools \
    --recreate_bed \
    --transcript_tsv /path/to/transcript_tsv
```

###Mandatory Arguments
Argument | Description
---|---
vcf | Path to germline VCF file of sample. For instance the germline vcf output from PURPLE. Calls should be wrt v37.
sample_t_id | The tumor sample ID of the run. Used for names of output files.
sample_r_id | The ref sample ID of the run.
version | The version of PEACH.
outputdir | Directory to write the output to.
panel | Path to a JSON file that contains the variants and haplotypes to test on.
vcftools | Path to VCFtools >= 0.1.14 (to allow for VCF v4.2).

###Optional Arguments
Argument | Default | Description
---|---|---
recreate_bed | False | To filter the VCF to the genes of interest, we use a transcript file and vcftools to filter on bed. Use this argument to regenerate the bed-file. If not given, the cached bed-file is used. The path to the cached bed file is "{path/to/panel/json}.bed".
transcript_tsv | None | If the bed file should be recreated, then this argument is required. This file should be a tsv file that describes transcripts for genes wrt v37, including the genes in the panel JSON.

####Datastore file locations
Panel:
  * Smaller panel for DPYD with haplotypes and haplotypes restricted to those in SOC tests (`/data/common/dbs/peach/panelfiles/min_DPYD.json`).
  * Panel with common DPYD haplotypes and variants (`/data/common/dbs/peach/panelfiles/DPYD.json`).

Transcript tsv: `/data/common/dbs/peach/all_genes.37.tsv`
