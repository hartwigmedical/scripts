# HMF_PGx README

HMF_PGx is a pharmacogenomics tool developed for the [Hartwig Medical Foundation pipeline](https://github.com/hartwigmedical/pipeline5).
It imports curated variants from a JSON file and inspects their presence the germline VCF. It creates two output files:
* [sample]_genotype.txt; contains on each line a determined genotype of the sample for a specific gene
* [sample]_calls.txt; contains all the variants that were used for testing and their respective calls and filters.
 
## Installation
If you want to run the code, please generate a local Python 3 venv and install the requirements:

```bash
$ python3 -m venv [path/to/new/virtual/environment, for example: ./pgx]
$ source [path/to/new/venv, for example: ./pgx/bin/activate]
(pgx) $ pip install -r requirements.txt
```

## Usage
Remember to source the virtualenv before running `main.py`.

####General usage
```(pgx) $ python main.py input.vcf.gz [optional args]```

####Arguments
* `vcf`: (Required) Path to germline VCF file of sample.
* `sampleID`: (Required) The sample ID of the run.
* `version`: (Required) The version of the tool.
* `--panel`: A file that contains the variants to test on. Different options are available for a panel file:
    * Curated JSON file (for example: `data/panelfiles/DPYD.json`). If this file is used, a hardcoded exceptions file is used to refactor hg19/hg38 differences.
    * TSV file where each row is a variant and containing the following columns: `chrom, hg19_start, gene, rsid`.
    * If no panel is given, a gene panel is queried from the PharmGKB API.
* `--outputdir`: Directory to write the output to. If no outputdir given, haplotypes will be written to stdout.
* `--vcftools`: Path to VCFtools >= 0.1.14 (to allow for VCF v4.2).
* `--requery`: If no panel is given in the --panel argument, this switch can be used to requery PharmGKB, instead of using the PharmGKB cache.
* `--recreate_bed`: To filter the VCF to the genes of interest, we use a transcript file and vcftools to filter on bed. Use this argument to regenerate the bed-file. If not given, the cached bed-file is used.
* `--tempdir`: Optional temp directory where intermediary VCF files are stored. Default = repodir/data
* `--sourcedir`: If source files are not loaded from repodir/data, an alternative location can be given.

