# PEACH

**P**harmacogenomic **E**valuator **A**nd **C**aller of **H**aplotypes (PEACH) is a pharmacogenomics tool developed for 
the [Hartwig Medical Foundation pipeline](https://github.com/hartwigmedical/pipeline5). 
It imports haplotypes and related variants from a curated JSON file, reports the presence of these variants in a 
germline VCF, and infers the simplest combination of haplotypes that explains the presence of these variants. 

It creates two output files:
* A file that contains the determined genotype of the sample for each gene in the JSON, expressed in terms of haplotypes.
* A file that contains all of the variants from the JSON file and their respective calls and filters.

## Contents

* [Installation](#installation)
* [Arguments](#arguments)
  + [Mandatory Arguments](#mandatory-arguments)
  + [Optional Arguments](#optional-arguments)
* [Input](#input)
  + [VCF](#vcf)
  + [JSON](#json)
  + [Transcript TSV](#transcript-tsv)
  + [Datastore file locations](#datastore-file-locations)
* [Output](#output)
  + [Genotype TSV file](#genotype-tsv-file)
  + [Calls TSV file](#calls-tsv-file)
* [Algorithm](#algorithm)

## Installation
If you want to run PEACH, please generate a local Python 3 venv and install the requirements:

```bash
$ python3 -m venv [path/to/new/virtual/environment, for example: ./peach]
$ source [path/to/new/venv, for example: ./peach/bin/activate]
(peach) $ pip install -r requirements.txt
```

## Arguments
Remember to source the virtualenv before running `main.py`.

#### Example Usage
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

### Mandatory Arguments
Argument | Description
---|---
vcf | Path to germline VCF file of sample. For instance the germline VCF output from PURPLE. Calls should be wrt v37.
sample_t_id | The tumor sample ID of the run. Used for names of output files.
sample_r_id | The ref sample ID of the run.
version | The version of PEACH.
outputdir | Directory to write the output to.
panel | Path to a JSON file that contains the variants and haplotypes to test on.
vcftools | Path to [VCFtools](http://vcftools.sourceforge.net/) >= 0.1.14 (to allow for VCF v4.2).

### Optional Arguments
Argument | Default | Description
---|---|---
recreate_bed | N/A | To filter the VCF to the genes of interest, we use a transcript file and VCFTools to filter on a bed file. Use this argument to regenerate the bed file. If not given, the cached bed-file is used. The path to the cached bed file is "{path/to/panel/json}.bed".
transcript_tsv | None | If the bed file should be recreated, then this argument is required. This file should be a tsv file that describes transcripts for genes wrt v37, including the genes in the panel JSON.

## Input
### VCF
PEACH has been designed to work with VCF files that follow the VCF Version 4.2 format, see 
[specification](https://samtools.github.io/hts-specs/VCFv4.2.pdf). In addition to the required fields, 
each row should contain an annotation field "ANN", as described in 
[annotation format specification](http://snpeff.sourceforge.net/VCFannotationformat_v1.0.pdf), that contains
the subsections "Gene Name" and "HGVS.c". One of the samples in the VCF should have a label 
equal to the `sample_r_id` argument,
and the "GT" sub-field for this sample should be included and filled in with diploid calls.

### JSON
For an example of a valid panel JSON (with fake data), see 
[example](https://github.com/hartwigmedical/scripts/blob/master/peach/src/test_resources/test_panel.json).
All fields in the example JSON are required. Additional fields are ignored. 
Relevant differences between the v37 and v38 reference sequences for a gene can be included as an entry in the "variants" field
of that gene where the "referenceAlleleV37" and "referenceAlleleV38" fields are different. The set of rs id's with such entries 
should be equal to the set of rs id's with entries in the "refSeqDifferenceAnnotations" field of that gene.

PEACH does not (properly) support panel JSON files that contain (partially) overlapping genes.
Variants in a panel JSON file are not allowed to (partially) overlap.

### Transcript TSV
TODO: write or give link

#### Datastore file locations 
(Only relevant for internal use)

Panel:
* Smaller panel for DPYD with haplotypes and haplotypes restricted to those in SOC tests (`/data/common/dbs/peach/panelfiles/min_DPYD.json`).
* Panel with common DPYD haplotypes and variants (`/data/common/dbs/peach/panelfiles/DPYD.json`).

Transcript tsv: `/data/common/dbs/peach/all_genes.37.tsv`

## Output
PEACH outputs two TSV files. One contains genotypes/haplotypes for each gene, the other contains calls for all of the variants from the panel JSON.
#### Genotype TSV file
Name: `[sample_t_id].peach.genotype.tsv`

Column | Example Value | Description
---|---|---
gene | DPYD | Gene for which this haplotype is called.
haplotype | *1_HOM | Haplotype from JSON, including whether it is homozygous (HOM) or heterozygous (HET). If no haplotype could be called, has value "Unresolved Haplotype".
function | No function | Functionality of this haplotype. Wild type has function "Normal Function". If no haplotype could be called, has value "Unknown Function".
linked_drugs | 5-Fluoracil;Capecitabine | Drugs for which this haplotype is relevant, separated by ";".
url_prescription_info | https://www.some_url.com/5-Fluoracil;https://www.some_other_url.com/Capecitabine | For each listed drug, a url with information on how to translate abnormal haplotype function into an appropriate treatment adjustement.
panel_version | DPYDpanel_v1.3 | Name and version of panel JSON. Both are taken from fields in the JSON.
repo_version | 1.0 | Version of PEACH.

#### Calls TSV file
Name: `[sample_t_id].peach.calls.tsv`

Column | Example Value | Description
---|---|---
gene | DPYD | Gene to which the variant is related.
chromosome | 1 | Chromosome of variant.
position_v37 | 98348885 | Position on chromosome wrt v37 reference genome.
position_v38 | 97883329 | Position on chromosome wrt v38 reference genome. If v37 info could not be translated into its v38 equivalent, has value "UNKNOWN".
ref_v37 | G | Reference allele wrt v37. 
ref_v38 | A | Reference allele wrt v38. If v37 info could not be translated into its v38 equivalent, has value "UNKNOWN".
allele1 | A | First of the called alleles. Order of alleles is lexicographical order.
allele2 | A | Second of the called alleles. Order of alleles is lexicographical order.
rsid | rs1801265 | Rs id(s) of variant. If more than one, then they are separated by ";". Taken from VCF if available. If not, taken from matching variant in panel JSON, if possible. If not, has value ".".
variant_annotation_v37 | 85C>T | Variant annotation wrt v37. See TODO for details.
filter_v37 | PASS | Has value PASS or NO_CALL. See TODO for details.
variant_annotation_v38 | REF_CALL | Variant annotation wrt v38. See TODO for details.
filter_v38 | NO_CALL | Has value PASS, NO_CALL, UNKNOWN or INFERRED_PASS. See TODO for details.
panel_version | DPYDpanel_v1.3 | Name and version of panel JSON. Both are taken from fields in the JSON.
repo_version | 1.0 | Version of PEACH.

TODO: better description special/missing values

TODO: Describe output VCF. Also in other parts of Readme?

## Algorithm
TODO: describe algorithm steps and provide details for each step

### Preparation
First, the panel JSON is loaded and checked for consistency. 

If the `--recreate_bed` argument was passed, 
then the bed file corresponding to the panel JSON is (re)created. 
To this end, for each gene in the panel JSON, corresponding positions are extracted from the transcript TSV such that 
the range between those start and end positions covers the entire gene.

### Get Variant Calls V37
Using VCFtools, the input VCF is filtered on the ranges in the bed file and on the sample name `sample_r_id`. 
The filtered VCF is read, and it is compared to the variants in the panel JSON file. 
Calls are ignored when none of the following are true:
* At least one of the rs id's of the call matches an rs id from the panel JSON.
* At least one of the covered positions of the call matches at least one of the covered positions of the variants in the panel JSON.
In this comparison, the *covered positions* of a call or a variant are the positions of the bases in the reference allele.

Let's call the remaining variants the *observed V37 calls*. 

For each variant in the panel JSON for which there are no matching calls in the observed V37 call, 
a call is created that is homozygously the reference allele, an *inferred V37 call*. 
More specifically, a call is added for a panel variant when there are no observed V37 calls such that either:
* The variant rs id matches one of the calls rs id's.
* The covered positions of the variant (partially) match the covered positions of the call.

The observed and unobserved V37 calls together form the list of calls that will be considered by PEACH.

### Add V38 Information to Calls
TODO: write




### Infer Haplotypes
TODO: write

### Produce Output
TODO: write

### Restrictions
TODO: write

PEACH does not support calling for multiple (partially) overlapping genes.
If one wishes to attain results for (partially) overlapping genes anyway,
split them across separate panel JSON files and run PEACH multiple times.

Variants in a panel JSON file are not allowed to (partially) overlap.

Differences in reference sequence between v37 and v38 should not be entered as MNV's, but as separate SNV's.