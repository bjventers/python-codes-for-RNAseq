'''This script parses junction file from Tophat and construct

splice junction objects, which are stored in dictionary using
a coordinate as a key.

Author: Likit Preeyanon
Email: preeyano@msu.edu

'''

import sys
import csv


class Junction(object):
    def __init__(self, chrom, start, end, coverage, strand, name):
        self.chrom = chrom
        self.start = start
        self.end = end
        self.coverage = coverage
        self.name = name
        self.strand = strand

    def __str__(self):
        return '%s, %s' % (self.getCoord(), self.name)

    def getCoord(self):
        return '%s:%d-%d' % (self.chrom, self.start, self.end)


def getJunction(row):


    chrom = row[0]
    name = row[3]
    coverage = int(row[4])
    strand = row[5]
    chromStart = int(row[1])
    blockStarts = [int(b) for b in row[-1].split(',')]
    blockSizes = [int(s) for s in row[-2].split(',')]

    for i in range(len(blockStarts) - 1):
        juncStart = blockStarts[i] + blockSizes[i] + chromStart
        juncEnd = blockStarts[i + 1] + chromStart
        junction = Junction(chrom, juncStart, juncEnd,
                                coverage, strand, name)
        yield junction


def parseJunctions(fileName):

    container = {}

    with open(fileName) as junctionFile:
        reader = csv.reader(junctionFile, dialect='excel-tab')
        try:
            for rowNum, row in enumerate(reader, start=1):
                if rowNum % 1000 == 0:
                    print('... %d' % rowNum)

                assert len(row) == 12, \
                '''A junction file from Topphat must
                
                contain exactly 12 columns

                '''

                blockCount = int(row[9])
                coverage = int(row[4])

                junctionNumber = 0
                for junction in getJunction(row):
                    container[junction.getCoord()] = junction
                    junctionNumber += 1

                assert junctionNumber == blockCount - 1, \
                    '''A number of junctions is less than a number of

                    block count by 1
                    
                    '''

        except csv.Error, e:
            sys.exit('file %s, line %d: %s' % (fileName, reader.line_num, e))

    return container

def findMatch(container1, container2):
    '''this function finds junctions that are common in two

    datasets.

    '''

    common = {}
    diff = {}
    for key in container1.keys():
        try:
            junction = container2[key]
        except KeyError:
            junc1 = container1[key]
            if junc1.coverage >= 10:
                diff[key] = junc1
                print(key, junc1.coverage)
        else:
            junc1 = container1[key]
            junc2 = container2[key]
            if junc1.coverage >= 10 and junc2.coverage >= 10:
                common[key] = junc1

    return common, diff


if __name__ == '__main__':

    introns = parseJunctions(sys.argv[1])
    print('\n')
    junctions1 = parseJunctions(sys.argv[2])
    print('\n')
    junctions2 = parseJunctions(sys.argv[3])
    print('\n')

    matchedJunctions, diffJunctions = findMatch(junctions1, junctions2)
    matchedModels, diffModels = findMatch(diffJunctions, introns)

    print('Dataset1 vs dataset 2')
    print(len(matchedModels), ' matched')
    print(len(diffModels), ' not matched')
    for key in matchedModels.keys():
        junc = matchedModels[key]
        print >> sys.stderr, key, junc.name, junc.coverage
