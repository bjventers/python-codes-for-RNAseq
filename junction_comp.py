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
        self.name = [name]
        self.strand = strand

    def __str__(self):
        return '%s, %s' % (self.getCoord(), self.name)

    def getCoord(self):
        return '%s:%d-%d' % (self.chrom, self.start, self.end)

class Model(object):
    def __init__(self, chrom, start, end, name, strand,
                            blockSizes, blockStarts):

        self.chrom = chrom
        self.start = int(start)
        self.end = int(end)
        self.blockSizes = [int(s) for s in blockSizes.split(',')]
        self.name = name
        self.strand = strand
        self.blockStarts = [int(s) for s in blockStarts.split(',')]

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
                    print >> sys.stderr, '... {0}'.format(rowNum)

                assert len(row) == 12, \
                '''A junction file from Topphat must
                
                contain exactly 12 columns

                '''

                blockCount = int(row[9])

                junctionNumber = 0
                for junction in getJunction(row):
                    try:
                        existingJunction = container[junction.getCoord()]
                    except KeyError:
                        container[junction.getCoord()] = junction
                    else:
                        existingJunction.name += junction.name
                            
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
        else:
            junc1 = container1[key]
            junc2 = container2[key]
            if junc1.coverage >= 10 and junc2.coverage >= 10:
                common[key] = junc1

    return common, diff

def scanJunctions(model, junctions1, junctions2, container):

    for i in range(len(model.blockStarts)):
        start = model.blockStarts[i] + model.start
        end = model.blockSizes[i] + start

        #key = '%s:%d-%d' % (model.chrom, start, end)
        juncKey = '%s:%d' % (model.chrom, end)

        try:
            junc1 = junctions1[juncKey]
        except KeyError:
            pass
        else:
            try:
                junc2 = junctions2[juncKey]
            except KeyError:
                pass
            else:
                junc1Ends = set(junc1)
                junc2Ends = set(junc2)
                diff = junc1Ends.difference(junc2Ends)
                if list(diff):
                    container[juncKey] = list(diff)

    return container


def buildJunctionDict(junctions):
    container = {}
    for k,v in junctions.iteritems():
        key = '%s:%d' % (v.chrom, v.start)
        try:
            junction = container[key]
        except KeyError:
            container[key] = [v.end]
        else:
            if v.end not in container[key]:
                container[key].append(v.end)

    return container


def findAlternativeSplicing(modelsFileName, junctions1, junctions2):

    container = {}

    with open(modelsFileName) as modelsFile:
        reader = csv.reader(modelsFile, dialect='excel-tab')
        try:
            for rowNum, row in enumerate(reader, start=1):
                if rowNum % 1000 == 0:
                    print >> sys.stderr, '... {0}'.format(rowNum)

                assert len(row) == 12, \
                '''A junction file from Topphat must
                
                contain exactly 12 columns

                '''

                chrom = row[0]
                start = row[1]
                end = row[2]
                name = row[3]
                strand = row[5]
                blockSizes = row[-2]
                blockStarts = row[-1]
                model = Model(chrom, start, end, name, strand,
                                    blockSizes, blockStarts)
                scanJunctions(model, junctions1, junctions2, container)

        except csv.Error, e:
            sys.exit('file %s, line %d: %s' % (fileName, reader.line_num, e))

    return container

if __name__ == '__main__':

    print >> sys.stderr, 'parsing models...'
    modelsFileName = sys.argv[1]
    models = parseJunctions(modelsFileName)
    print >> sys.stderr, '\n'

    print >> sys.stderr, 'parsing junctions...'
    junctions1 = buildJunctionDict(parseJunctions(sys.argv[2]))
    print >> sys.stderr, '\n'

    print >> sys.stderr, 'parsing junctions...'
    junctions2 = buildJunctionDict(parseJunctions(sys.argv[3]))
    print >> sys.stderr, '\n'

    print >> sys.stderr, 'searching for an alternative splicing...'
    altSplicing = findAlternativeSplicing(modelsFileName, junctions1, junctions2)

    for k in altSplicing:
        for j in altSplicing[k]:
            key = '%s-%d' % (k, j)
            try:
                transcript = models[key]
            except KeyError:
                pass
            else:
                print('{0}\t{1}'.format(key, ','.join(transcript.name)))
