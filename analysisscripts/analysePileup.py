
PATH = "/Users/peterpriestley/hmf/sliceCPCT02010267/CPCT02010267.pileup"
PILEUP_FILENAME = ""
sorted_coordinates = [(17,55999856),(17,60000147)]


def readPileup(aPath,aFileName,aPileupCoordinates):
    aPileupCoordinates.reverse()
    myCoordinate = aPileupCoordinates.pop()
    filteredPileup = []
    with open(aPath+aFileName, 'r') as f:
        for line in f:
            line = line.strip('\n')
            a = [x for x in line.split('\t')]
            while aPileupCoordinates and (int(a[0])>myCoordinate[0] or (int(a[0])==myCoordinate[0] and int(a[1]) > myCoordinate[1])):
                myCoordinate = aPileupCoordinates.pop()
            if int(a[0]) == myCoordinate[0] and int(a[1]) == myCoordinate[1]:
                filteredPileup.append(a)
    return filteredPileup

pileupOutput = readPileup(PATH,PILEUP_FILENAME,sorted_coordinates)
print pileupOutput

