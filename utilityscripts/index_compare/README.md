# index_compare

Code to check all combinations of any number of barcode index collections. Determines the i7 and i5
distance for every possible combination. Result for each is FAIL in case total distance is
lower than the allowed mismatches during conversion * 2 + 1 and otherwise OK.

### Example TSV input dual index collection
    #i7name i7seq   i5name  i5seq
    bc1i7   CTGATCG bc1i5   GCGCATAT
    bc2i7   CTGATCG bc2i5   GCGCATAA

See also [the standard dual index collection](hmf_indexes_dual_8.tsv)

### Example output
    #result dist1 dist2 idx1_tag1 idx1_seq1 idx1_tag2 idx1_seq2 idx2_tag1 idx2_seq1 idx2_tag2 idx2_seq2
    OK 5 3 F498_LCH588 CTGCGGAT CZ-0379-LL_A AGATAACC IDT8_i7_384 CGACCATT IDT8_I5_384 TGATAGGC