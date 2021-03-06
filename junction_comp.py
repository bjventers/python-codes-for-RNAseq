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
        self.chromStart = int(start)
        self.chromEnd = int(end)
        self.blockSizes = [int(s) for s in blockSizes.split(',')]
        self.name = name
        self.strand = strand
        self.blockStarts = [int(s) for s in blockStarts.split(',')]

    def __str__(self):
        return '%s, %s' % (self.getCoord(), self.name)

    def getCoord(self):
        return '%s:%d-%d' % (self.chrom, self.start, self.end)


class Exon(object):
    def __init__(self, chrom, start, end,
                    prevStart, prevEnd,
                    nextStart, nextEnd,
                    geneName):

        self.chrom = chrom
        self.start = start
        self.end = end

        if prevStart and prevEnd:
            prevExon = '%s:%d-%d' % (self.chrom, prevStart, prevEnd)
            self.prevExons = set([prevExon])
        else:
            self.prevExons = set([])

        if nextStart and nextEnd:
            nextExon = '%s:%d-%d' % (self.chrom, nextStart, nextEnd)
            self.nextExons = set([nextExon])
        else:
            self.nextExons = set([])

        self.geneName = geneName


class ModelJunction(object):
    def __init__(self, chrom, start, end, event):
        self.chrom = chrom
        self.start = start
        self.end = end
        self.event = event

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

        juncStart = '%s:%d' % (model.chrom, end)

        try:
            junc1 = junctions1[juncStart]
        except KeyError:
            pass
        else:
            try:
                junc2 = junctions2[juncStart]
            except KeyError:
                pass
            else:
                junc1Ends = set(junc1)
                junc2Ends = set(junc2)
                diff = junc1Ends.difference(junc2Ends)
                if list(diff):
                    container[juncStart] = list(diff)

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


def findAlternativeSplicing(junctions1, junctions2):

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
            sys.exit('file %s, line %d: %s' % (modelsFileName, reader.line_num, e))

    return container


def parseGeneModels(modelsFileName):

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
                yield model

        except csv.Error, e:
            sys.exit('file %s, line %d: %s' % (modelsFileName,
                                                reader.line_num, e))

def groupExons(genes, exons, isoform):
    geneName, isonum = isoform.name.split('.')

    for i in range(len(isoform.blockStarts)):
        start = isoform.blockStarts[i] + isoform.chromStart
        end = isoform.blockSizes[i] + start

        if i == 0:
            prevStart = None
            prevEnd = None
        else:
            prevStart = isoform.blockStarts[i-1] + isoform.chromStart
            prevEnd = prevStart + isoform.blockSizes[i-1]
            prevExon = '%s:%d-%d' % (isoform.chrom, prevStart, prevEnd)

        if i == len(isoform.blockStarts)-1:
            nextStart = None
            nextEnd = None
        else:
            nextStart = isoform.blockStarts[i+1] + isoform.chromStart
            nextEnd = nextStart + isoform.blockSizes[i+1]
            nextExon = '%s:%d-%d' % (isoform.chrom, nextStart, nextEnd)

        exonCoord = '%s:%d-%d' % (isoform.chrom, start, end)
        
        try:
            e = exons[exonCoord]
        except KeyError:
            exon = Exon(isoform.chrom, start, end,
                        prevStart, prevEnd,
                        nextStart, nextEnd,
                        geneName)
            exons[exonCoord] = exon
        else:
            if prevStart: e.prevExons.add(prevExon)
            if nextStart: e.nextExons.add(nextExon)

        try:
            gene = genes[geneName]
        except KeyError:
            genes[geneName] = set([exonCoord])
        else:
            if exonCoord not in genes[geneName]:
                genes[geneName].add(exonCoord)


def identifyJunctions(genes, exons):
    junctions = {}
    def getStart(coord):
        chrom, startEnd = coord.split(':')
        start, end = [int(j) for j in startEnd.split('-')]
        return start

    for k in genes.keys():
        allExons = sorted(genes[k], key=lambda x: getStart(x))
        for i in range(len(allExons)):
            exonCoord = allExons[i]
            exon = exons[exonCoord]

            '''
            if len(exon.nextExons) > 1:
                for nextExonCoord in exon.nextExons:
                    nextExon = exons[nextExonCoord]

                    if len(exons[nextExonCoord].prevExons) > 1:
                        altEvent = 'skippedExon'
                    else:
                        altEvent = 'alternativeSpliceSite'

                    junction = ModelJunction(exon.chrom, exon.end,
                                        nextExon.start, altEvent)

                    juncCoord = '%s:%d-%d' % (exon.chrom,
                                                exon.end,
                                                nextExon.start)
                    junctions[juncCoord] = junction
            '''

            if len(exon.prevExons) > 1:
                for prevExonCoord in exon.prevExons:
                    prevExon = exons[prevExonCoord]

                    if len(prevExon.nextExons) > 1:
                        for e in prevExon.nextExons:
                            ex = exons[e]
                            if ex.end < exon.start:
                                altEvent = 'skippedExon'
                                break
                    else:
                        if prevExon.prevExons == set([]):
                            altEvent = 'alternativeSplicing'
                        else:
                            for e in exon.prevExons:
                                otherPrevExon = exons[e]
                                if otherPrevExon.end != prevExon.end and \
                                    otherPrevExon.end > prevExon.start:
                                    altEvent = 'alternativeSpliceSite'
                                    break
                                else:
                                    altEvent = None

                    junction = ModelJunction(exon.chrom, exon.end,
                                        prevExon.start, altEvent)

                    juncCoord = '%s:%d-%d' % (exon.chrom,
                                                prevExon.end,
                                                exon.start)
                    junctions[juncCoord] = junction

            if len(exon.prevExons) == 1:
                for prevExonCoord in exon.prevExons:
                    prevExon = exons[prevExonCoord]

                    if len(prevExon.nextExons) > 1:
                        for e in prevExon.nextExons:
                            ex = exons[e]
                            print >> sys.stderr, ex.start, ex.end
                            if ex.start != exon.start:
                                if ex.nextExons == set([]) or \
                                    exon.nextExons == set([]):
                                    altEvent = 'alternativeSplicing'
                                else:
                                    altEvent = 'alternativeSpliceSite'
                                break
                            else:
                                altEvent = None
                    else:
                        altEvent = None

                    junction = ModelJunction(exon.chrom, prevExon.end,
                                        exon.start, altEvent)
                    juncCoord = '%s:%d-%d' % (exon.chrom,
                                                prevExon.end,
                                                exon.start,
                                                )
                    try:
                        j = junctions[juncCoord]
                    except KeyError:
                        junctions[juncCoord] = junction

            if len(exon.nextExons) == 1:
                altEvent = None

                for nextExonCoord in exon.nextExons:
                    nextExon = exons[nextExonCoord]
                    junction = ModelJunction(exon.chrom, exon.end,
                                        nextExon.start, altEvent)

                    juncCoord = '%s:%d-%d' % (exon.chrom,
                                                exon.end,
                                                nextExon.start)
                    try:
                        j = junctions[juncCoord]
                    except KeyError:
                        junctions[juncCoord] = junction

    #print >> sys.stderr, junctions.keys()
    return junctions

if __name__ == '__main__':

    modelsFileName = sys.argv[1]

    genes = {}
    exons = {}

    print >> sys.stderr, 'Tagging model junctions ...'
    for model in parseGeneModels(modelsFileName):
        groupExons(genes, exons, model)

    print >> sys.stderr, '%d genes, %d exons' % (len(genes), len(exons))

    modelJunctions = identifyJunctions(genes, exons)

    print >> sys.stderr, '%d junctions' % len(modelJunctions)

    skipped = [k for k in modelJunctions.keys() if modelJunctions[k].event == \
                'skippedExon']

    altss = [k for k in modelJunctions.keys() if modelJunctions[k].event == \
                'alternativeSpliceSite']

    noAlt = [k for k in modelJunctions.keys() if not modelJunctions[k].event]


    print >> sys.stderr, '\n'
    print >> sys.stderr, 'parsing junctions 1 ...'
    junctions1 = parseJunctions(sys.argv[2])
    print >> sys.stderr, '\n'

    print >> sys.stderr, 'parsing junctions 2 ...'
    junctions2 = parseJunctions(sys.argv[3])
    print >> sys.stderr, '\n'

    diff1 = set(junctions1.keys()).difference(set(junctions2.keys()))
    diff2 = set(junctions2.keys()).difference(set(junctions1.keys()))

    print >> sys.stderr, '%d skipped exons' % len(skipped)
    print >> sys.stderr, '%d skipped exons' % len(altss)
    print >> sys.stderr, '%d skipped exons' % len(noAlt)
    print >> sys.stderr, 'Junctions1 diff junctions2 %d junctions' % len(diff1)
    print >> sys.stderr, 'Junctions2 diff junctions1 %d junctions' % len(diff2)

    for k in diff1:
        try:
            j = modelJunctions[k]
        except KeyError:
            pass
        else:
            print('{0}\t{1}\t1\t{2}'.format(k, j.event, junctions1[k].coverage))

    print >> sys.stderr, k
    for k in diff2:
        try:
            j = modelJunctions[k]
        except KeyError:
            pass
        else:
            print('{0}\t{1}\t2\t{2}'.format(k, j.event, junctions2[k].coverage))
