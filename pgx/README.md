# HMF_PGx README

HMF_PGx is a pharmacogenomics tool developed for the [Hartwig Medical Foundation pipeline](https://github.com/hartwigmedical/pipeline5).
It imports haplotypes and related variants from a curated JSON file, and inspects the presence of these variants in the germline VCF. 

It creates two output files:
* [sample]_genotype.txt; contains on each line a determined genotype of the sample for a specific gene.
* [sample]_calls.txt; contains all the variants from the JSON file and their respective calls and filters.
 
## Installation
If you want to run the code, please generate a local Python 3 venv and install the requirements:

```bash
$ python3 -m venv [path/to/new/virtual/environment, for example: ./pgx]
$ source [path/to/new/venv, for example: ./pgx/bin/activate]
(pgx) $ pip install -r requirements.txt
```

## Usage
Remember to source the virtualenv before running `main.py`.

####Example usage
```
(pgx) $ python main.py input.vcf.gz \
    COLO829T \
    COLO829R \
    1.0 \
    /path/to/outputdir/ \
    /path/to/panel.json \
    /path/to/vcftools \
    --recreate_bed \
    --sourcedir /path/to/sourcedir
```

####Arguments
* `vcf`: (Required) Path to germline VCF file of sample. For instance the germline vcf output from PURPLE. Calls should be wrt v37.
* `sample_t_id`: (Required) The tumor sample ID of the run.
* `sample_r_id`: (Required) The ref sample ID of the run.
* `version`: (Required) The version of the tool.
* `outputdir`: (Required) Directory to write the output to.
* `panel`: (Required) A Curated JSON file that contains the variants and haplotypes to test on. Different options are available for a panel file:
    * Panel with common DPYD haplotypes and variants (`/data/panelfiles/DPYD.json`).
    * Smaller panel for DPYD with haplotypes and haplotypes restricted to those in SOC tests (`/data/panelfiles/min_DPYD.json`).
* `vcftools`: (Required) Path to VCFtools >= 0.1.14 (to allow for VCF v4.2).
* `--recreate_bed`: (Optional, default=False) To filter the VCF to the genes of interest, we use a transcript file and vcftools to filter on bed. 
  Use this argument to regenerate the bed-file. If not given, the cached bed-file is used. 
  The path to the cached bed file is the path to the panel file, except with ".json" replaced by ".bed"
* `--transcript_tsv`: (Optional, default=/data/common/dbs/pgx/all_genes.37.tsv) 
  If the bed file should be recreated, then this argument is required. 
  This file should be a tsv file that describes transcripts for genes wrt v37, including the genes in the panel JSON.



