# PEACH

**P**harmacogenomic **E**valuator **A**nd **C**aller of **H**aplotypes (PEACH) is a pharmacogenomics tool developed for 
the [Hartwig Medical Foundation pipeline](https://github.com/hartwigmedical/pipeline5). 
It imports haplotypes and related variants from a curated JSON file, reports the presence of these variants in a 
germline VCF, and infers the simplest combination of haplotypes that explains the presence of these variants. 

The two main output files are:
* A file that contains the determined genotype of the sample for each gene in the JSON, expressed in terms of haplotypes.
* A file that contains calls for all positions of variants in the JSON file, including annotation and filters wrt both v37 and v38 reference genomes.

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
  + [Preparation](#preparation)
  + [Get Variant Calls V37](#get-variant-calls-v37)
  + [Annotate Calls with Panel Information](#annotate-calls-with-panel-information)
  + [Infer Haplotypes](#infer-haplotypes)
  + [Examples](#examples)
  + [Restrictions](#restrictions)

## Installation
PEACH has been designed to work with Python 3.6.

If you want to run PEACH, please generate a local Python 3.6 venv and install the requirements:

```bash
$ python3.6 -m venv [path/to/new/virtual/environment, for example: ./peach]
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

The calls in the VCF should be with respect to a v37 reference genome.

### JSON
For an example of a valid panel JSON (with fake data), see 
[example](https://github.com/hartwigmedical/scripts/blob/master/peach/src/test_resources/test_panel.json).
All fields in the example JSON are required. Additional fields are ignored. 
Relevant differences between the v37 and v38 reference sequences for a gene should be included as an entry in the "variants" field
of that gene where the "referenceAlleleV37" and "referenceAlleleV38" fields are different. The set of rs id's with such entries 
has to be equal to the set of rs id's with entries in the "refSeqDifferenceAnnotations" field of that gene.
This ensures that the variant annotation for variants at these locations can be translated to v38 properly

PEACH does not (properly) support panel JSON files that contain (partially) overlapping genes.
Variants in a panel JSON file are not allowed to (partially) overlap.

### Transcript TSV
TODO: write or provide link

#### Datastore file locations 
(Only relevant for internal use)

Panel:
* Smaller panel for DPYD with haplotypes and haplotypes restricted to those in SOC tests (`/data/common/dbs/peach/panelfiles/min_DPYD.json`).
* Panel with common DPYD haplotypes and variants (`/data/common/dbs/peach/panelfiles/DPYD.json`).

Transcript tsv: `/data/common/dbs/peach/all_genes.37.tsv`

## Output
PEACH outputs two TSV files. One contains genotypes/haplotypes for each gene, the other contains calls for all variants from the panel JSON.
PEACH also outputs a VCF file that is a filtered version of the input VCF.
### Genotype TSV file
Name: `[sample_t_id].peach.genotype.tsv`

Column | Example Value | Description
---|---|---
gene | DPYD | Gene for which this haplotype is called.
haplotype | *1_HOM | Haplotype from JSON, including whether it is homozygous (HOM) or heterozygous (HET). If no haplotype could be called, has value "Unresolved Haplotype".
function | No function | Functionality of this haplotype. Wild type has function "Normal Function". If no haplotype could be called, has value "Unknown Function".
linked_drugs | 5-Fluoracil;Capecitabine | Drugs for which this haplotype is relevant, separated by ";".
url_prescription_info | https://www.some_url.com/5-Fluoracil;https://www.some_other_url.com/Capecitabine | For each listed drug, a url with information on how to translate abnormal haplotype function into an appropriate treatment adjustment. Separated by ";".
panel_version | DPYDpanel_v1.3 | Name and version of panel JSON. Both are taken from fields in the JSON.
repo_version | 1.0 | Version of PEACH.

### Calls TSV file
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
rsid | rs1801265 | Rs id(s) of variant. If more than one, then they are separated by ";". Taken from VCF if available. If not, taken from matching variant in panel JSON, if match exists. If not, has value ".".
variant_annotation_v37 | 85C>T | Variant annotation wrt v37. See [Get Variant Calls V37](#get-variant-calls-v37) for details.
filter_v37 | PASS | Has value PASS or NO_CALL. See [Get Variant Calls V37](#get-variant-calls-v37) for details.
variant_annotation_v38 | REF_CALL | Variant annotation wrt v38. See [Annotate Calls with Panel Information](#annotate-calls-with-panel-information) for details.
filter_v38 | NO_CALL | Has value PASS, NO_CALL, UNKNOWN or INFERRED_PASS. See [Annotate Calls with Panel Information](#annotate-calls-with-panel-information) for details.
panel_version | DPYDpanel_v1.3 | Name and version of panel JSON. Both are taken from fields in the JSON.
repo_version | 1.0 | Version of PEACH.

### Filtered VCF
Name: `[sample_t_id].filtered.vcf`

The result of filtering the input VFC by sample id `sample_r_id` and by a bed file with ranges for genes from the JSON.
See [Get Variant Calls V37](#get-variant-calls-v37) for a more detailed explanation.

## Algorithm
Haplotypes are commonly defined wrt a v38 reference genome. 
Since PEACH accepts VCF files wrt v37 as input, this requires a translation of v37 calls to v38 calls.
PEACH extracts the required knowledge of the differences between these reference genomes from the information in the panel JSON.

In broad strokes, PEACH does the following:
* Extract relevant calls wrt from VCF, where relevance is determined by the panel JSON. 
  These calls will be with respect to a v37 reference genome.
* Use information from the panel JSON to translate this set of calls wrt v37 into the corresponding set of calls wrt v38, where possible.
* For each gene:
  + Determine for each variant how often each alt allele occurs.
  + Determine the unique simplest combination of haplotypes that completely explains that combination of alt alleles and counts.
If there is no unique simplest combination of haplotypes that completely explains the combination of alt alleles and counts, then declare "Unresolved Haplotype".
* Create output files.

### Preparation
First, the panel JSON is loaded and checked for consistency. 

If the `--recreate_bed` argument was passed, 
then the bed file corresponding to the panel JSON is (re)created. 
To this end, for each gene in the panel JSON, corresponding positions are extracted from the transcript TSV such that 
the range between those start and end positions covers the entire gene.

### Get Variant Calls V37
Using VCFtools, the input VCF is filtered on the ranges in the bed file and on the sample name `sample_r_id`. 
The resulting file is included in the output as `[sample_t_id].peach.filtered.vcf`.

The filtered VCF is read, and it is compared to the variants in the panel JSON file. 
Calls are ignored when none of the following are true:
* At least one of the rs id's of the call matches an rs id from the panel JSON.
* At least one of the covered positions of the call matches at least one of the covered positions of the variants in the panel JSON.
In this comparison, the *covered positions* of a call or variant are the positions of the bases in the reference allele.

Let's call the remaining variants the *observed v37 calls*. 

For each variant in the panel JSON for which there are no matching calls in the observed v37 call, 
a call is created that is homozygously the reference allele, an *inferred v37 call*. 
More specifically, a call is added for a panel variant when there are no observed v37 calls such that either:
* The variant rs id matches one of the calls rs id's.
* The covered positions of the variant (partially) match the covered positions of the call.

The observed and inferred v37 calls together form the list of calls that will be considered by PEACH, 
which we will call the *combined v37 calls*.

The annotation and filter of the combined v37 calls are determined in the following way:

Condition | Variant Annotation V37 | Filter V37
---|---|---
Homozygous alt or heterozygous observed call. | From HGVS.c field in ANN in VCF | PASS  
Homozygous reference observed call. | REF_CALL | PASS
Inferred call. | REF_CALL | NO_CALL

### Annotate Calls with Panel Information
For each of the combined v37 calls, an attempt is made to find a variant in the panel JSON that has the same v37 position and reference allele.

If successful:
* If the rs id of the call has not been set, then it is set to the value from the matching variant.
* The reference allele and position wrt v38 are determined from the matching variant.

If unsuccessful:
* Set reference allele and position wrt v38 to "UNKNOWN".

Also, the correct annotation and filter wrt v38 are determined according to the following table, 
where an asterisk (*) indicates that any value is allowed:

Matches variant in panel JSON | Type of call | Call is homozygous reference wrt v38 | Reference alleles v37 and v38 are identical | All alleles are reference with v37 or v38 | Variant Annotation V38 | Filter V38
---|---|---|---|---|---|---
False|*|*|*|*|Variant Annotation V37 + "?"|UNKNOWN
True|Inferred|True|*|*|REF_CALL|NO_CALL
True|Inferred|False|*|*|From "refSeqDifferenceAnnotations" field in panel JSON|INFERRED_PASS
True|Observed|True|*|*|REF_CALL|PASS
True|Observed|False|True|*|Variant Annotation V37|PASS
True|Observed|False|False|True|From "refSeqDifferenceAnnotations" field in panel JSON|PASS
True|Observed|False|False|False|Variant Annotation V37 + "?"|UNKNOWN

Note that an asterisk does not indicate that every value is actually possible. 
For instance, calls that do not match any variants from the panel JSON can only be observed calls, not inferred calls.
Furthermore, if there is no matching variant in the panel JSON, then PEACH does not know what the correct v38 reference allele is, 
so it would be unknown whether the reference alleles wrt v37 and v38 are identical. 

Let's call these calls with both v37 and v38 details *full calls*.

### Infer Haplotypes
The goal is to find the simplest combination of haplotypes that explains the called variants. 

Sometimes, more than one combination of haplotypes could explain these calls. 
As an example, consider the DPYD gene and two variants for that gene: c.1905+1G>A and c.1627A>G. 
Separately, these variants form the haplotypes *2A and *5, and the haplotype *2B consists of both of these variants together.
A haplotype can contain multiple variants if these variants have a tendency to be inherited together.
If each variant is called once and if all of these variants and haplotypes are included in the panel JSON, then, to take this combined inheritance into account, 
PEACH prefers to call the haplotype combination as *2B_HET/*1_HET and not as *2A_HET/*5_HET. 
If you want PEACH to call *2A_HET/*5_HET instead of *2B_HET/*1_HET in this situation, then simply don't include *2B in the panel JSON.

To make this more precise, define the *length* of a haplotype combination as the total number of non-wild-type
haplotypes in the combination, where homozygous haplotype calls are counted as 2. 
PEACH will always attempt to call the unique haplotype combination of minimum length that explain all of the variant calls.
If there are no haplotype combinations that explain all of the variant calls, 
or if there is more than one combination of the same minimum length, 
then the haplotype combination for that gene is called as "Unresolved Haplotype".

The only valid haplotype combination of length 0 is the homozygous wild type haplotype,
valid haplotype combinations of length 1 always include precisely one heterozygous wild type haplotype call,
and valid haplotype combinations of length at least 2 do not contain any calls for the wild type haplotype.

#### Haplotype Calling Algorithm 
Haplotypes are called for each gene separately. First, collect the full calls that correspond to that gene. 
Then, extract the alt alleles wrt v38 from these full calls, 
and count the number of times each combination of position (v38) and alt allele (v38) occurs.
Use recursive descent to determine all haplotype combinations that perfectly explain all of these variants.
If there are no such combinations, then no haplotype combination can be called for this gene. 
If such combinations do exist, then the next step is to determine the length of
each valid haplotype combination. Find the minimum length of the valid haplotype combinations, 
and select the haplotype combinations whose length is equal to this minimum. 
If precisely one such haplotype combination exists, then this combination will be called for this gene.
If more than one haplotype combination of minimum length exists, then no haplotype combination is called for this gene.

### Examples
The data in these examples will be the completely fictional.
The examples will focus on fairly "standard" situations, and they will exclude all information that is not necessary to understand these situations. 
For details on non-standard situations, see the more detailed subsections of the [Algorithm](#algorithm) section.

Suppose that the panel JSON contains the following variants and haplotypes for the fictional gene FAKE,
and that FAKE is the only gene in the panel JSON.

Rs Id|Reference Allele V37|Reference Allele V38|V38 Annotation for Reference Sequence Difference
---|---|---|---
rs1|A|A|N/A
rs2|TA|GC|c.6543GC>TA
rs3|GG|GG|N/A

Haplotype|Variants (Rs Id: Variant Allele wrt v38)
---|---
*1 (wild type)|None
*2|rs1: T
*3|rs2: TA
*4|rs3: G
*5|rs1: T, rs3: G

#### No Calls
If there are no calls wrt v37 in the VCF, then the full calls are:

Rs Id|Allele1|Allele2|Variant Annotation V37|Filter V37|Variant Annotation V38|Filter V38
---|---|---|---|---|---|---
rs1|A|A|REF_CALL|NO_CALL|REF_CALL|NO_CALL
rs2|TA|TA|REF_CALL|NO_CALL|c.6543GC>TA|INFERRED_CALL
rs3|GG|GG|REF_CALL|NO_CALL|REF_CALL|NO_CALL

The only valid haplotype combination that explains these variants is *3_HOM, 
so this is the haplotype combination that is called for FAKE.

#### Homozygous Wild Type
Suppose that the v37 calls are the following:

Rs Id|Allele1|Allele2|Variant Annotation V37|Filter V37
---|---|---|---|---
rs1|A|A|c.8483A>T|PASS
rs2|GC|GC|c.6543TA>GC|PASS
rs3|GG|GG|c.4838GG>G|PASS

In this case, the full calls are:

Rs Id|Allele1|Allele2|Variant Annotation V37|Filter V37|Variant Annotation V38|Filter V38
---|---|---|---|---|---|---
rs1|A|A|REF_CALL|PASS|REF_CALL|PASS
rs2|GC|GC|c.6543TA>GC|PASS|REF_CALL|PASS
rs3|GG|GG|REF_CALL|PASS|REF_CALL|PASS

The only valid haplotype combination is *1_HOM, so this haplotype combination is called for FAKE.

#### Heterozygous Wild Type
Suppose that the v37 calls are the following:

Rs Id|Allele1|Allele2|Variant Annotation V37|Filter V37
---|---|---|---|---
rs1|A|A|c.8483A>T|PASS
rs2|TA|GC|c.6543TA>GC|PASS
rs3|GG|GG|c.4838GG>G|PASS

The full calls are:

Rs Id|Allele1|Allele2|Variant Annotation V37|Filter V37|Variant Annotation V38|Filter V38
---|---|---|---|---|---|---
rs1|A|A|REF_CALL|PASS|REF_CALL|PASS
rs2|GC|TA|c.6543TA>GC|PASS|c.6543GC>TA|PASS
rs3|GG|GG|REF_CALL|PASS|REF_CALL|PASS

The only valid haplotype combination for FAKE is *3_HET/*1_HET, so this haplotype combination is called.

#### Multiple Valid Haplotypes
Suppose that the v37 calls are the following:

Rs Id|Allele1|Allele2|Variant Annotation V37|Filter V37
---|---|---|---|---
rs1|A|T|c.8483A>T|PASS
rs2|GC|GC|c.6543TA>GC|PASS
rs3|G|G|c.4838GG>G|PASS

The resulting full calls are:

Rs Id|Allele1|Allele2|Variant Annotation V37|Filter V37|Variant Annotation V38|Filter V38
---|---|---|---|---|---|---
rs1|A|T|c.8483A>T|PASS|c.8483A>T|PASS
rs2|GC|GC|c.6543TA>GC|PASS|REF_CALL|PASS
rs3|G|G|c.4838GG>G|PASS|c.4838GG>G|PASS

The valid haplotype combinations are *2_HET/*4_HOM and *4_HET/*5_HET. 
These combinations have lengths 3 and 2, respectively, so the second combination is preferred.
The called haplotype combination for FAKE is *4_HET/*5_HET.

### Restrictions
PEACH does not support calling for multiple (partially) overlapping genes.
If one wishes to attain results for (partially) overlapping genes anyway,
split them across separate panel JSON files and run PEACH multiple times.

Variants in a panel JSON file are not allowed to (partially) overlap.