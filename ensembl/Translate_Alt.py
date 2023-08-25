

#!/usr/bin/python

TranscriptFile="inputFiles/MajorityVoteAtLeast2_subversion_adapt.txt"
CurrentEnsemblFile="inputFiles/ensembl_trans_exon_data.csv"
refSeqFile="inputFiles/refseqtableSubversionComp2votes_noAlt.csv"
EnsemblGenesFile="inputFiles/ensembl_gene_data.csv"

EnsemblGenesfh=open(EnsemblGenesFile,'r')
CanonicalTransScriptIds=dict()
GeneNameToEnsemblGeneId=dict()
EnsemblGeneDict=dict()
EnsemblGeneIdToGeneName=dict()
refSeqDict = dict()
MajorityVoteDict = dict()

transcriptCount=987654321


def printRefSeq(geneName,transcriptCount):

  arr=refSeqDict[geneName]

  refSeqId = arr[0]
  refSeqTranscriptId = arr[1]
  refSeqTranscriptId = refSeqTranscriptId.split(".")[0]
  if (refSeqTranscriptId == "NM_021059") or (refSeqTranscriptId == "NM_022148"): #No overlap with gene start to gene end (H3C14 and CRLF2)
    return 0
  if not (refSeqTranscriptId in MajorityVoteDict):
    return 0
  refSeqchr = arr[2]
  refSeqStrand = arr[3]
  if refSeqStrand == "-":
    refSeqStrand="-1"
  else:
    refSeqStrand="1"
  refSeqTransStart = arr[4]
  refSeqTransEnd = arr[5]
  refSeqcdsStart = arr[6]
  refSeqcdsEnd = arr[7]
  refSeqExonCount = arr[8]
  refSeqExonStarts = arr[9]
  refSeqExonStarts=refSeqExonStarts.split(",")

  refSeqExonEnds = arr[10]
  refSeqExonEnds = refSeqExonEnds.split(",")
  refSeqScore = arr[11]
  refSeqName2 = arr[12]
  refSeqStartStat = arr[13]
  refSeqEndStat = arr[14]
  refSeqExonFrames = arr[15]
  refSeqExonFrames = refSeqExonFrames.split(",")

#  print(arr)
#  transcriptCount+=1
  if refSeqName2  in EnsemblGeneDict:
    arr = EnsemblGeneDict[refSeqName2]
    for i in range(int(refSeqExonCount)):
      if i == int(refSeqExonCount)-1:
        refSeqExonEndPhase = "-1"
      else:
        refSeqExonEndPhase = str(int(refSeqExonFrames[i+1]))
      if refSeqStrand=="-1":
        print(arr[0]+","+arr[1]+","+ refSeqStrand+","+str(transcriptCount)+","+refSeqTranscriptId+","+"protein_coding"+","+str(int(refSeqTransStart)+1)+","+refSeqTransEnd+","+str(int(refSeqExonCount)-i)+","+str(int(refSeqExonStarts[i])+1)+","+refSeqExonEnds[i]+","+refSeqExonFrames[i]+","+refSeqExonEndPhase+","+str(int(refSeqcdsStart)+1)+","+refSeqcdsEnd)
      else:
         print(arr[0]+","+arr[1]+","+ refSeqStrand+","+str(transcriptCount)+","+refSeqTranscriptId+","+"protein_coding"+","+str(int(refSeqTransStart)+1)+","+refSeqTransEnd+","+str(i+1)+","+str(int(refSeqExonStarts[i])+1)+","+refSeqExonEnds[i]+","+refSeqExonFrames[i]+","+refSeqExonEndPhase+","+str(int(refSeqcdsStart)+1)+","+refSeqcdsEnd)
  else:
    print("ERROR: " + refSeqName2 +" NOT FOUND")


first=True

MajorityVotefh = open(TranscriptFile,'r')

for lines in MajorityVotefh:
  lines=lines.rstrip("\n")
  arr = lines.split(".")
  MajorityVoteDict[arr[0]]=1


refSeqfh = open(refSeqFile,'r')


first = True

for lines in refSeqfh:
  if first:
   first= False
   continue

  lines=lines.rstrip('\n')
  arr=lines.split("\t")
  refSeqDict[arr[12]]=arr
first = True

for lines in EnsemblGenesfh:

  if not first:
    lines=lines.rstrip('\n')
    arr=lines.split(",")
    GeneNameToEnsemblGeneId[arr[1]]=arr[0]
    EnsemblGeneIdToGeneName[arr[0]]=arr[1]
  else:
    first = False


CurrentEnsemblfh = open(CurrentEnsemblFile,'r')

first=True
prev='Empty'

for lines in CurrentEnsemblfh:
  lines=lines.rstrip("\n")
  if not first:
    arr=lines.split(",")
    EnsemblGeneDict[EnsemblGeneIdToGeneName[arr[0]]]=arr
    if arr[0] != prev and prev !='Empty' and EnsemblGeneIdToGeneName[prev] in refSeqDict:
       printRefSeq(EnsemblGeneIdToGeneName[prev],transcriptCount)
       transcriptCount+=1
       prev=arr[0]
    elif arr[0]!= prev:
      prev = arr[0]
  print(lines)
  first = False
#print(EnsemblGeneDict)
#Transcriptfh=open(TranscriptFile,'r')


#for lines in Transcriptfh:
#  lines=lines.rstrip("\n")
#  arr=lines.split("\t")
#  Gene = arr[0]
#  majorityVote = arr[9]
#  hartwigTranscript = arr[12]
#  print(Gene,majorityVote,hartwigTranscript)





