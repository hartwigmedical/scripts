#!/usr/local/bin/python

chromosomeLength = {}
chromosomeLength['1'] = 249250622
chromosomeLength['10'] = 135534749
chromosomeLength['11'] = 135006518
chromosomeLength['12'] = 133851897
chromosomeLength['13'] = 115169880
chromosomeLength['14'] = 107349541
chromosomeLength['15'] = 102531394
chromosomeLength['16'] = 90354755
chromosomeLength['17'] = 81195212
chromosomeLength['18'] = 78077250
chromosomeLength['19'] = 59128985
chromosomeLength['2'] = 243199374
chromosomeLength['20'] = 63025522
chromosomeLength['21'] = 48129897
chromosomeLength['22'] = 51304568
chromosomeLength['3'] = 198022431
chromosomeLength['4'] = 191154277
chromosomeLength['5'] = 180915261
chromosomeLength['6'] = 171115068
chromosomeLength['7'] = 159138664
chromosomeLength['8'] = 146364023
chromosomeLength['9'] = 141213432
chromosomeLength['MT'] = 16571
chromosomeLength['X'] = 155270561
chromosomeLength['Y'] = 59373567

def intChrom(chrom):
    if chrom == 'X':
        return 23
    elif chrom == 'Y':
        return 24
    elif chrom == 'MT':
        return 25
    elif str(chrom)[:2] == 'GL':
        return 100
    elif str(chrom)[:2] == 'NC':
        return 100
    else:
        return int(chrom)