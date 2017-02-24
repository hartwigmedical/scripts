#!/usr/bin/env python
import heapq
import os
import argparse

class vcfMerge():

    def __init__(self,outputFile, minCountThreshold):

        try:
            self._heap = []
            self._output_file = open(outputFile, 'w+')
            self.minCountThreshold = minCountThreshold

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
        return "\t"

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
            read_line_split = read_line.split("\t")
            heapq.heappush(self._heap, (str(self.posToNumber(read_line_split))+"\t"+read_line_split[4], file__))


    def writeToOutput(self,positionAlt,count):
        if count >= self.minCountThreshold:
            self._output_file.write(self.numberToPos(float(positionAlt.split("\t")[0])) + self._delimiter_value() \
                                    + positionAlt.split("\t")[1]+ self._delimiter_value() + str(count) + "\n")

def getVCFList(path,suffixMatch):
    vcfList = []
    for path, subdirs, files in os.walk(path):
        for name in files:
            if name[-len(suffixMatch):] == suffixMatch:
                vcfList.append(os.path.join(path, name))
    return vcfList

if __name__ == '__main__':
    parser = argparse.ArgumentParser(formatter_class=lambda prog: argparse.HelpFormatter(prog, max_help_position=100, width=200))
    required_named = parser.add_argument_group('required named arguments')
    required_named.add_argument('-s', '--suffixMatch', help="file suffix to match eg: 'xxx.vcf'", required=True)
    required_named.add_argument('-m', '--minCountThreshold', help='minCount to add to PON output.  eg: 2', required=True)
    required_named.add_argument('-p', '--path', help='directory to search for matching files', required=True)
    required_named.add_argument('-o', '--outputFile', help='output file name', required=True)
    args = parser.parse_args()

    files = getVCFList(args.path,args.suffixMatch)
    merger = vcfMerge(args.outputFile,int(args.minCountThreshold))
    merger.merge(files)



