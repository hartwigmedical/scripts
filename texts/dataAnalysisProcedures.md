# Data Analysis Procedures

### Patient report procedure

For tumor biopsies within the CPCT and DRUP studies, Hartwig Medical Foundation provides centers with a "patient report" that contains the most relevant somatic findings. The reported variants are based on the final list of somatic variants from our analysis pipeline but in addition have to pass certain "rules" before they end up in the report:

- Consensus rule: only variants with enough support from the various somatic callers are reported (see details below)
- Consequence rule: only variants located in the canonical transcript region of certain genes that also have a certain predicted effect are reported (see details below)

**Details Consensus rule:**  
- include SNV variant if found by 3 or more callers (out of 4 total)
- include INDEL variant if found by 2 or more callers (out of 3 total)
- include variant if located in the "high confidence region" based AND is called by at least 2 callers AND is found in COSMIC or not in DBSNP
- include variant if located within "CPCT slicing panel" genes (75 MBase, 65 genes)

**Details Consequence rule:**  
- the list of 118 genes (covering about 15.5 MBases) can be found in any patient report
- only variants with one of the following consequences/effects are repoted: transcript ablation, splice_acceptor_variant
splice_donor_variant, stop_gained, frameshift_variant, stop_lost, initiator_codon_variant, start_lost, transcript_amplification, inframe_insertion, inframe_deletion, missense_variant, splice_region_variant, incomplete_terminal_codon_variant
- we use the "canonical transcript" definition of ensembl: the canonical transcript for a gene is set according to the following hierarchy: 1. Longest CCDS translation with no stop codons. 2. If no (1), choose the longest Ensembl/Havana merged translation with no stop codons. 3. If no (2), choose the longest translation with no stop codons. 4. If no translation, choose the longest non-protein-coding transcript.

**Mutational load:**  
The total number of missense variants is reported as the "mutational load". The definition of missense is as follows: a sequence variant, that changes one or more bases, resulting in a different amino acid sequence but where the length is preserved. For all consequence definitions see the [calculated variant consequences table](http://www.ensembl.org/info/genome/variation/predicted_data.html) at Ensembl.

**Links:**  
[Source code Consensus Rule Filter](https://github.com/hartwigmedical/hmftools/tree/master/consensus-rule-filter)  
[Source code Patient Reporter](https://github.com/hartwigmedical/hmftools/tree/master/patient-reporter)

-----

