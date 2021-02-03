select * from somaticVariant
where (
(chromosome = 7 and gene = "MET" and position < 116411903 and position > 116411888 and reported =0)
or (chromosome = 7 and gene = "MET" and position < 116412058 and position > 116412043 and reported =0)

or (chromosome = 17 and gene = "TP53" and position < 7571720 and position > 7571713 and reported =0)
or (chromosome = 17 and gene ="TP53" and position > 7573008 and position < 7573015 and reported =0)

or (chromosome = 17 and gene = "TP53" and position < 7573927 and position > 7573920 and reported =0)
or (chromosome = 17 and gene ="TP53" and position > 7574033 and position < 7574040 and reported =0)

or (chromosome = 17 and gene = "TP53" and position < 7576853 and position > 7576846 and reported =0)
or (chromosome = 17 and gene ="TP53" and position > 7576926 and position < 7576933 and reported =0)

or (chromosome = 17 and gene = "TP53" and position < 7577019 and position > 7577012 and reported =0)
or (chromosome = 17 and gene ="TP53" and position > 7577155 and position < 7577162 and reported =0)

or (chromosome = 17 and gene = "TP53" and position < 7577499 and position > 7577492 and reported =0)
or (chromosome = 17 and gene ="TP53" and position > 7577608 and position < 7577615 and reported =0)

or (chromosome = 17 and gene = "TP53" and position < 7578177 and position > 7578170 and reported =0)
or (chromosome = 17 and gene ="TP53" and position > 7578289 and position < 7578296 and reported =0)

or (chromosome = 17 and gene = "TP53" and position < 7578371 and position > 7578364 and reported =0)
or (chromosome = 17 and gene ="TP53" and position > 7578554 and position < 7578561 and reported =0)

or (chromosome = 17 and gene = "TP53" and position < 7579312 and position > 7579305 and reported =0)
or (chromosome = 17 and gene ="TP53" and position > 7579590 and position < 7579597 and reported =0)

or (chromosome = 17 and gene = "TP53" and position < 7579700 and position > 7579693 and reported =0)
or (chromosome = 17 and gene ="TP53" and position > 7579721 and position < 7579728 and reported =0))
and sampleId in ("XXX");