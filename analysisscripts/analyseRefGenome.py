#!/usr/local/bin/python
import random

randomFastaLength = 40000000
refgenome = "/Users/peterpriestley/hmf/data/refgenomes/chr21/chr21.fasta"
randomPath = "/Users/peterpriestley/"
randomFileName = "random.fasta"

chars = {}
chars[1] = 'A'
chars[2] = 'T'
chars[3] = 'G'
chars[4] = 'C'

def countKmers(aRefGenome,aKmers,aKmerLength):
    print 'reading:',aRefGenome,"for kmers length",aKmerLength
    with open(aRefGenome, 'r') as f:
        myKmer = ""
        while True:
            ch = f.read(1)
            if not ch: break
            if ch == '>':
                f.readline()
            elif ch != "\n" and ch != 'N':
                myKmer = myKmer + ch
                if len(myKmer) > aKmerLength:
                    myKmer = myKmer[1:]
                if len(myKmer) == aKmerLength:
                    if aKmers.has_key(myKmer):
                        aKmers[myKmer] += 1
                    else:
                        aKmers[myKmer] = 1

    return aKmers


# LOGIC FOR COUNTING REPEATS
# keep kmer of length = maxRepeatLength * 2
# For each base processed
    # For j = 1 to maxRepeatLength
        # repeatTracker(j) = [sequence of kmer len = j, # of repeats observed (0+), current pos in repeat sequence]
        # if base = next expected base in repeat
            # if current pos < j
                # then increment currentPos
            # else
                # increment # of repeats
        # else
            # pop sequence + append to sequence
            # reset # of repeats = 0
            # Write to repeate dictionary:   [len, # of repeats] += 1

def countRepeats(aRefGenome,aRepeats,aMaxRepeatLength):
    print 'reading:',aRefGenome,"for repeats"
    with open(aRefGenome, 'r') as f:
        myKmer = ""
        while True:
            ch = f.read(1)
            if not ch: break
            if ch == '>':
                f.readline()
            elif ch != 'N':
                pass  # NEED TO ADD LOGIC HERE
            elif ch != "\n" and ch != 'N':
                pass  # NEED TO ADD LOGIC HERE
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

if __name__ == "__main__":
    kmers = {}
    kmersRandom = {}
    kmerLength = 7
    kmers = countKmers(refgenome,kmers,kmerLength)
    for kmer, kmerCount in sorted(kmers.items()):
        print kmer, ":REAL:", kmerCount
    generateRandomFasta(randomPath,randomFileName,randomFastaLength)
    kmersRandom = countKmers(randomPath + randomFileName,kmersRandom,kmerLength)
    for kmer, kmerCount in sorted(kmersRandom.items()):
        print kmer, ":RAND:", kmerCount

    repeats = {}
    maxRepeatLength = 10
    repeats = countRepeats(refgenome,repeats,maxRepeatLength)
