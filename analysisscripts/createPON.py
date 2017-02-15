
import heapq
import os

#### TO DO ###########
# Min, Max and Median AF?
# Look for exact SNP MATCH instead of just pos???

class vcfMerge():

    minCountThreshold = 2

    def __init__(self,path):

        try:
            self._heap = []
            self._output_file = open(path+'PON.tsv', 'w+')

        except Exception, err_msg:
            print "Error while creating Merger: %s" % str(err_msg)

    def merge(self, input_files):
        try:
            open_files = []
            [open_files.append(open(file__, 'r')) for file__ in input_files]

            [self.readFirstVariant(file__) for file__ in open_files]

            previousCount = 0
            smallest = self._heap[0]

            while(self._heap):

                previous = smallest[0]
                smallest = heapq.heappop(self._heap)

                if previous < smallest[0]:
                    self.writeToOutput(previous,previousCount)
                    previousCount = 1
                else:
                    previousCount += 1

                read_line = smallest[1].readline()
                self.pushVariantToHeap(read_line, smallest[1])

            self.writeToOutput(smallest[0], previousCount)

            [file__.close() for file__ in open_files]
            self._output_file.close()

        except Exception, err_msg:
            print "Error while merging: %s" % str(err_msg)

    def _delimiter_value(self):
        return "\n"

    def posToNumber(self,split):
        if split[0] == 'X':
            chrom = 23
        elif split[0] == 'Y':
            chrom = 24
        elif split[0] == 'MT':
            chrom = 25
        else:
            chrom = int(split[0])
        return chrom * 1e9 + int(split[1])

    def numberToPos(self, number):
        if int(number/1e9) == 23:
            return "X\t"+ str(int(number%1e9))
        elif int(number/1e9) == 24:
            return "Y\t" + str(int(number % 1e9))
        elif int(number / 1e9) == 25:
            return "MT\t" + str(int(number % 1e9))
        else:
            return str(int(number/1e9)) + "\t" + str(int(number%1e9))

    def readFirstVariant(self,file__):
        read_line = file__.readline()
        while read_line[0] == "#":
            read_line = file__.readline()
        self.pushVariantToHeap(read_line,file__)

    def pushVariantToHeap(self,read_line,file__):
        if(len(read_line) != 0):
            heapq.heappush(self._heap, (self.posToNumber(read_line.split("\t")), file__))


    def writeToOutput(self,posNumber,count):
        if count >= self.minCountThreshold:
            self._output_file.write(self.numberToPos(posNumber) + "\t" + str(count) + self._delimiter_value())

def getVCFList(path):
    files = []
    for x in os.listdir(path):
        if x[-4:] == ".vcf":
            files.append(path + x)
    return files

if __name__ == '__main__':
    #path = "/Users/peterpriestley/hmf/analyses/PON/"
    path = "/Users/peterpriestley/hmf/analyses/170117_PMC_analysis/"
    files = getVCFList(path)
    merger = vcfMerge(path)
    merger.merge(files)



