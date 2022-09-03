# index_compare

Contains code to check all combinations of any number of barcode index collections. Determines the i7 and i5 distance for every possible combination. Result for each is FAIL in case total distance is
lower than the allowed mismatches during bcl2fastq conversion * 2 + 1 and otherwise OK.

The [whitelist file](whitelisted_dual_indexes.tsv) contains dual indexes of length 8 that are safe to use in combination with [our indexes](hmf_indexes_dual_8.tsv).

### Example TSV input dual index collection
    #i7name i7seq   i5name  i5seq
    bc1i7   CTGATCG bc1i5   GCGCATAT
    bc2i7   CTGATCG bc2i5   GCGCATAA

### Example output
    #result dist1 dist2 idx1_tag1 idx1_seq1 idx1_tag2 idx1_seq2 idx2_tag1 idx2_seq1 idx2_tag2 idx2_seq2
    OK 5 3 internal-1-i7 CGACCATT internal-1-i5 TGATAGGC external-1-i7 CTGCGGAT external-1-i5 AGATAACC
    FAIL 2 2 internal-2-i7 TTCCAAGG internal-i5 CTACAAGG external-2-i7 TACCGAGG external-2-i5 CAACAATG
 