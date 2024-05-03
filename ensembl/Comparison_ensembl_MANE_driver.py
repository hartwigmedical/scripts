
#!/usr/bin/python3



def readDriverGenePanel(fn):
  genesInDriverGenePanel=[]
  first=True
  with open(fn,'r') as fh:
    for lines in fh:
      lines=lines.rstrip('\n')
      if first == True:
        first = False
        continue
      arr = lines.split('\t')
      genesInDriverGenePanel.append(arr[0])
  return genesInDriverGenePanel

def getGeneIds(genesInDriverGenePanel,geneListFile):
  EnsemblToGene = dict()
  GeneToEnsembl = dict()  
  first = True
  with open(geneListFile,'r') as fh:
    for lines in fh:
      lines = lines.rstrip("\n")
      if first == True:
        first = False
        continue
      arr = lines.split(',')
      if arr[1] in genesInDriverGenePanel:
        EnsemblToGene[arr[0]]=arr[1]
        GeneToEnsembl[arr[1]]=arr[0]
  return(GeneToEnsembl)

def getCanonicalTranscripts(genesInDriverGenePanel,EnsemblToGene,transcriptFile):
  EnsemblToCanonical = dict()
  first = True
  with open(transcriptFile,'r') as fh:
    for lines in fh:
      lines= lines.rstrip('\n')
      if first == True:
        first = False
        continue
      arr = lines.split(',')
      if (arr[1] in genesInDriverGenePanel) and (arr[3]=='true'):
        EnsemblToCanonical[arr[0]] = arr[2]            
  return EnsemblToCanonical

def getRefSeqToGeneHartwig(file):
  geneToRefSeq = dict()
  first =True
  with open(file,'r') as fh:
    for lines in fh:
      if first == True:
        first = False
        continue
      lines=lines.rstrip('\n')
      arr = lines.split('\t')
      geneToRefSeq[arr[0]] = arr[1]
  return geneToRefSeq

genesInDriverGenePanel = readDriverGenePanel('DriverGenePanel.38.tsv')


canonicalTranscripts = getCanonicalTranscripts(genesInDriverGenePanel,0,'38/ensembl_trans_amino_acids.csv')
canonicalTranscripts37 = getCanonicalTranscripts(genesInDriverGenePanel,0,'37/ensembl_trans_amino_acids.csv')

GeneToEnsembl = getGeneIds(genesInDriverGenePanel,'38/ensembl_gene_data.csv')
GeneToEnsembl37 = getGeneIds(genesInDriverGenePanel,'37/ensembl_gene_data.csv')
geneToRefSeq = getRefSeqToGeneHartwig('hartwigRefSeqList.txt')

maneFile ='MANE.GRCh38.v1.2.summary.txt'
with open(maneFile,'r') as fh:
  status = ''
  for lines in fh:
    lines = lines.rstrip('\n')
    arr = lines.split('\t')
    currHg38='NA'
    currHg37='NA'
    if arr[3] in genesInDriverGenePanel:
      if GeneToEnsembl[arr[3]] in canonicalTranscripts:
        currHg38 = canonicalTranscripts[GeneToEnsembl[arr[3]]]
      if arr[3] in genesInDriverGenePanel:
        currHg37 = canonicalTranscripts37[GeneToEnsembl37[arr[3]]]
      if arr[3] in geneToRefSeq:
        if arr[5] == geneToRefSeq[arr[3]]:
          status = 'MATCH'
        else:
          status = 'MISMATCH'
        print(arr[3]+'\t'+arr[5]+'\t'+arr[7] + '\t'+ geneToRefSeq[arr[3]] +'\t'+currHg37 + '\t' + currHg38+ '\t' +arr[9]+ '\t' + status)
      else:
        print(arr[3]+'\t'+arr[5]+'\t'+arr[7] + '\t'+ 'NA' + '\t'+currHg37+'\t'+currHg38+ '\t' +arr[9]+ '\t' + 'NA')
     
