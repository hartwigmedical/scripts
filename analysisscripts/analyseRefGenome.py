#!/usr/local/bin/python
import random
import chromosomeDefinition

randomFastaLength = 40000000
refgenome = "/Users/peterpriestley/hmf/data/refgenomes/chr21/chr21.fasta"
refgenome = "/Users/peterpriestley/hmf/data/refgenomes/Homo_sapiens.GRCh37.GATK.illumina/Homo_sapiens.GRCh37.GATK.illumina.fasta"


#refgenome = "/Users/peterpriestley/hmf/testGenome"

randomPath = "/Users/peterpriestley/"
randomFileName = "random.fasta"

def factors(n):
    return set(reduce(list.__add__,
                ([i, n//i] for i in range(1, int(n**0.5) + 1) if n % i == 0)))

chars = {}
chars[1] = 'A'
chars[2] = 'T'
chars[3] = 'G'
chars[4] = 'C'

def countKmers(aRefGenome,aKmers,aKmerLength):
    print 'reading:',aRefGenome,"kmers of length",aKmerLength
    with open(aRefGenome, 'r') as f:
        myKmer = ""
        while True:
            ch = f.read(1)
            if not ch: break
            if ch == '>':
                f.readline()
                myKmer = ""
            elif ch != 'N':
                myKmer = ""
            elif ch != "\n":
                myKmer = myKmer + ch
                if len(myKmer) > aKmerLength:
                    myKmer = myKmer[1:]
                if len(myKmer) == aKmerLength:
                    if aKmer.has_key(myKmer):
                        aKmer[myKmer] += 1
                    else:
                        akmer[myKmer] = 1
    return aKmers

def countRepeats(aRefGenome,aRepeats,aMaxRepeatLength):
    print 'reading:', aRefGenome, "repeats up to length", aMaxRepeatLength
    repeatTracker = []  #[[repeatSequence,numRepeats,currentPosInRepeat]]
    for i in range(maxRepeatLength):
        repeatTracker.append (["", 1, 0])
        i += 1
    with open(aRefGenome, 'r') as f:
        while True:
            ch = f.read(1)
            if not ch: break
            elif ch == 'N' or ch == '>':
                for i in range(maxRepeatLength):
                    repeatTracker[i] = ["", 1, 0]
                    i += 1
                for i in range(maxRepeatLength):
                    if repeatTracker[i][1] > 1:
                        if aRepeats.has_key((i + 1, repeatTracker[i][1])):  # add to aRepeats if numRepeats > 1
                            aRepeats[(i + 1, repeatTracker[i][1])] += 1
                        else:
                            aRepeats[(i + 1, repeatTracker[i][1])] = 1
                if ch == '>': # new chromosome.  reset
                    f.readline()
            elif ch != "\n":
                for i in range(maxRepeatLength):
                    if len(repeatTracker[i][0]) < i + 1:
                        repeatTracker[i][0]= repeatTracker[i][0] + ch
                    elif repeatTracker[i][0][(repeatTracker[i][2])%(i+1)] == ch: #potential repeat

                        if repeatTracker[i][2] == i:  # complete repeat - increment numRepeats
                            repeatTracker[i][1] += 1
                            repeatTracker[i][2] = (repeatTracker[i][2]+1)%(i+1)
                        else:  # increment current position
                            repeatTracker[i][2] += 1
                    else:  # no repeat - pop  sequence and log
                        if repeatTracker[i][1] > 1:
                            if aRepeats.has_key((i + 1, repeatTracker[i][1])):  # add to aRepeats if numRepeats > 1
                                aRepeats[(i + 1, repeatTracker[i][1])] += 1
                            else:
                                aRepeats[(i + 1, repeatTracker[i][1])] = 1
                        repeatTracker[i][0] = repeatTracker[i][0][repeatTracker[i][2]+1:] + repeatTracker[i][0][:repeatTracker[i][2]] + ch
                        repeatTracker[i][1] = 1
                        repeatTracker[i][2] = 0
    return aRepeats

def removeSubRepeats(aRepeats):

    for myRepeat,myRepeatCount in sorted(aRepeats.items()):
        myfactors = factors(myRepeat[0])
        for myfactor in myfactors:
            if myfactor > 0 and myfactor < myRepeat[0]:
                for i in range(myRepeat[0]*myRepeat[1]/myfactor,myRepeat[0]*myRepeat[1]/myfactor+myfactor):
                    if aRepeats.has_key((myfactor,i)):
                        aRepeats[myRepeat] -= aRepeats[(myfactor,i)]
                        print myRepeat[0],myRepeat[1],myRepeatCount,myfactor,i
    return aRepeats


def generateRandomFasta(aPath, aFilename,aLength):
    print 'Creating random Fasta file'
    with open(aPath + aFilename,'w') as f:
        f.write(">myFakeFasta\n")
        for i in range(aLength):
            a = random.randrange(1, 5)
            f.write(chars[a])
            if i % 60 == 59:
                f.write("\n")
    print 'File created:',aPath+aFilename

def GCContent(aRefGenome,windowSize):
    print 'reading:',aRefGenome,"trinucleotide frequencies"
    trinucleotides = {}
    myPos = 0
    with open(aRefGenome, 'r') as f:
        myKmer = ""
        while True:
            ch = f.read(1)
            if not ch: break
            if ch == '>':
                myChrom = (f.readline()[:-1])
                myPos = 0
            elif ch != "\n":
                myKmer = myKmer + ch
                if len(myKmer) > windowSize:
                    myKmer = ch
                    myPos = myPos + windowSize
                if len(myKmer) == windowSize:
                    if myKmer.count('N') < windowSize:
                        myGCContent =float((myKmer.count('G')+ myKmer.count('C')))/(windowSize - myKmer.count('N'))
                    else:
                        myGCContent = -1
                    print myChrom, myPos, myGCContent



def trinucleotideFrequencies(aRefGenome):
    print 'reading:',aRefGenome,"trinucleotide frequencies"
    trinucleotides = {}
    with open(aRefGenome, 'r') as f:
        myKmer = ""
        while True:
            ch = f.read(1)
            if not ch: break
            if ch == '>':
                f.readline()
                myKmer = ""
            elif ch == 'N':
                myKmer = ""
            elif ch != "\n":
                myKmer = myKmer + ch
                if len(myKmer) > 3:
                    myKmer = myKmer[1:]
                if len(myKmer) == 3:
                    if myKmer[1]=="C" or myKmer[1]=="T":
                        if trinucleotides.has_key(myKmer):  # add to aRepeats if numRepeats > 1
                            trinucleotides[myKmer] += 1
                        else:
                            trinucleotides[myKmer] = 1
    print len(trinucleotides)
    for trinucleotide,count in sorted(trinucleotides.items()):
        print trinucleotide,":",count

GCContent(refgenome,2000000)
#trinucleotideFrequencies(refgenome)


if __name__ == "__main__":
    kmers = {}
    kmersRandom = {}
    kmerLength = 7
    #kmers = countKmers(refgenome,kmers,kmerLength)
    #for kmer, kmerCount in sorted(kmers.items()):
    #    print kmer, ":REAL:", kmerCount
    #generateRandomFasta(randomPath,randomFileName,randomFastaLength)
    #kmersRandom = countKmers(randomPath + randomFileName,kmersRandom,kmerLength)
    #for kmer, kmerCount in sorted(kmersRandom.items()):
    #    print kmer, ":RAND:", kmerCount

    #repeats = {}   # {(kmerLength,repeatLength):count}
    #maxRepeatLength = 100
    #repeats = countRepeats(refgenome,repeats,maxRepeatLength)

    #for repeat, repeatCount in sorted(repeats.items()):
    #    print repeat[0],":",repeat[1],":",repeatCount
    #repeats = removeSubRepeats(repeats)
    #for repeat, repeatCount in sorted(repeats.items()):
    #    print repeat[0], ":", repeat[1], ":", repeatCount