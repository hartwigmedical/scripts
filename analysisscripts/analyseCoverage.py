
#!/usr/local/bin/python
import numpy as np
import matplotlib.pyplot as plt


sample = ('Benchmark',"/Users/peterpriestley/hmf/depthCoverage/","CPCT11111111T_dedup_WGSMetrics.txt")
sample2 = ('Weird',"/Users/peterpriestley/hmf/depthCoverage/","CPCT02020248T_dedup_WGSMetrics.txt")




def loadCoverageStats(aSample):
    print "reading coverage file. . .\n"
    mySampleName,myPath,myFilename = aSample
    coverage = []
    beginData = False
    with open(myPath + myFilename, 'r') as f:
        #i=0
        for line in f:
            line = line.strip('\n')
            a = [x for x in line.split('\t')]
            if beginData == True and a[0] != '':
                coverage.append((int(a[0]),int(a[1])))
            if a[0] == 'coverage':
                beginData = True

    return coverage


def printStatistics(aCoverageStats,aSample,aFigure):

    mySampleName,myPath,myFilename = aSample
    sumCoverage = 0.0
    for i in aCoverageStats:
        sumCoverage += i[1]
    print sumCoverage

    xValues = tuple(x[0] for x in aCoverageStats if (int(x[0])>-1 and int(x[0])<251))
    yValues = tuple((y[1]/sumCoverage) for y in aCoverageStats if (int(y[0])>-1 and int(y[0])<251))
    cumSumValues = 1-np.cumsum(yValues)


    fig = plt.figure(aFigure)
    ax = fig.add_subplot(111)
    rects1 = ax.bar(xValues, cumSumValues, color='r')

    ax.set_ylabel('% of Bases')
    ax.set_xlabel('Depth coverage')
    ax.set_title('Picard Depth Coverage: ' + mySampleName)

    plt.show()

if __name__ == "__main__":
    coverageStats = []
    coverageStats = loadCoverageStats(sample)
    printStatistics(coverageStats,sample,1)
    coverageStats = loadCoverageStats(sample2)
    printStatistics(coverageStats, sample2,2)
