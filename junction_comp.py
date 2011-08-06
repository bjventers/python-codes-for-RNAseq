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


def getJunctionStartEnd(row):
    chromStart = int(row[1])
    blockStarts = [int(b) for b in row[-1].split(',')]
    blockSizes = [int(s) for s in row[-2].split(',')]
    juncStart = blockStarts[0] + blockSizes[0] + chromStart
    juncEnd = blockStarts[-1] + chromStart
    return juncStart, juncEnd


def parseJunctions(fileName):

    junctions = {}

    with open(fileName) as junctionFile:
        reader = csv.reader(junctionFile, dialect='excel-tab')
        try:
            for row in reader:
                assert len(row) == 12, \
                '''
                    A junction file from Topphat must contain
                    exactly 12 columns
                '''

                chrom = row[0]
                name = row[3]
                coverage = int(row[4])
                strand = row[5]
                juncStart, juncEnd = getJunctionStartEnd(row)

                junction = Junction(chrom, juncStart, juncEnd,
                                        coverage, strand, name)

                junctions[junction.getCoord()] = junction

        except csv.Error, e:
            sys.exit('file %s, line %d: %s' % (fileName, reader.line_num, e))

    return junctions


if __name__ == '__main__':

    junctions = parseJunctions(sys.argv[1])
    for k, v in junctions.items():
        print(v)
