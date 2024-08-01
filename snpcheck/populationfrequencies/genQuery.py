#!/usr/bin/python

infile = 'SNPs.vcf'

bq=dict()

first=True

with open(infile,'r') as fh:
  
  for lines in fh:

    if first == True:
      first=False
      continue
    if lines[0] == '#':
      continue
    lines=lines.rstrip('\n')
    arr = lines.split('\t')
    chr = arr[0]
    if not (chr in bq):
      bq[chr] = "bq  query --max_rows 1500  --format json  --use_legacy_sql=false 'SELECT * FROM bigquery-public-data.gnomAD.v3_genomes__chr"
      bq[chr]+=str(chr)+ ' WHERE '
      bq[chr] += 'start_position = ' + str(int(arr[1])-1)+ ' OR '

#     print(lines)
    else:
      bq[chr] += 'start_position = ' + str(int(arr[1])-1)+ ' OR '
  for chr in bq:
    bq[chr] = bq[chr].rstrip(" OR ")
    bq[chr]+=" LIMIT 11000;' > bq_output_" + chr+".json"
    print(bq[chr])
