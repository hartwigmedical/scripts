import chromosomeDefinition as cd

def loadCNVFile(aBEDFile):
    myBed = []
    with open(aBEDFile, 'r') as f:
        chrom = "0"
        pos = 0
        f.readline()
        for line in f:
            line = line.strip('\n')
            splitLine = line.split(',')
            if (splitLine[2] > chrom or int(splitLine[3]) > pos):
                chrom = splitLine[2]
                pos = int(splitLine[3])
                print splitLine[2],chrom, splitLine[3], pos
                myBed.append(splitLine)
    return myBed



CSV_PATH="/Users/peterpriestley/hmf/analyses/CNVanalyses/"
CSV_FILENAME="copyNum.csv"
CNV = loadCNVFile(CSV_PATH + CSV_FILENAME)
