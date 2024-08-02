#!/usr/bin/python


#East asian
#Ashkenazi Jewish
#European (non-Finnish)
#Other
#Amish
#Latino/Admixed American
#Middle Eastern
#Afrian/African American
#XX
#YY

import json
import logging


POPS = [ "afr","ami","amr","asj","eas","fin","nfe","sas","oth","female","male"]
FIELDS = ['AN']
FIELDSab = ['AF']

header = 'chr\tstart_position\tref\talt\tAF\t'

for F in FIELDSab:
    for P in POPS:
        header+=F+"_"+P+"\t"
#header = header.rstrip('\t')

for F in FIELDS:
    for P in POPS:
        header+=F+"_"+P+"\t"
header = header.rstrip('\t')

currentSNPs = dict()

json_dict = dict()
file_probes="SNPs.vcf"


def genKey(chrJunran, posJunran):
    posJunran = str(int(posJunran)-1)
    return(chrJunran+"_"+posJunran)

def extractSNPs(filename,chrom):
    fh =open(filename)
    first = True
    with open('SNPs_addpopfreq.tsv','w') as output_fh:
        for lines in fh:
            lines=lines.rstrip("\n")

            if lines[0:2 ] ==  '##':
               output_fh.write(lines+'\n')
               continue
            if first:
                output_fh.write(lines+'\t'+header+"\n")
                first = False
            else:
                arr = lines.split("\t")
                if genKey(arr[0],arr[1]) in json_dict:
                    output_fh.write(lines+'\t'+json_dict[genKey(arr[0],arr[1])]+"\n")
                elif arr[0] == chrom:
                    logging.error((genKey(arr[0],arr[1]) + "not found"))

def extractjson(filename,chrom):
    fh = open(filename,'r')
    allVariants = json.load(fh)
    res=''
    for variant in allVariants:
        skip_variant = False
        ref = variant['reference_bases']
        chr = variant['reference_name']
        chr=chr[3:]
        alternate_bases=variant["alternate_bases"]
        start_position = variant["start_position"]
        for ab in alternate_bases:
            AF = ab['AF']
            if float(AF)<0.3 or float(AF)>0.7:
                 logging.warning("SKIPPING: " + start_position+ ", VAF: " + AF )
                 skip_variant = True
                 continue
            res = AF+'\t'

            vep = ab['vep']
            alt = ''
            for t in vep:
                if alt == '':
                    alt = (t['allele'])
                elif t['allele'] != alt:
                    logging.warning("different alt in transcript")
                    logging.warning(print(t))
            for F in FIELDSab:
                for P in POPS:
                    key = F+"_" +P
                    print(key)
                    res += str(ab[key]) + "\t"
        for F in FIELDS:
            for P in POPS:
                key = F+"_" +P
                res += variant[key] + "\t"
        res = res.rstrip("\t")
        print(json_dict)
        print(chr+"_"+start_position)
        if skip_variant == False and res != '':
            json_dict[chr+"_"+start_position] = chr+"\t"+start_position+'\t'+ref+'\t'+alt+'\t'+res
    return 0

for i in range(1,19):
  if (i == 3) or (i==4):
    continue
  chrom = str(i)
  
  x=extractjson("bq_output_" +chrom+".json",chrom)
#chrom ='X'
#x=extractjson("bq_output_" +chrom+".json",chrom)
print(json_dict)
extractSNPs(file_probes,chrom)
