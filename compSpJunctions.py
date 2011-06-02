import psl_parser
import sys

def construct(tName, tStarts, blockSizes, junctions):
    h = 0
    while h < len(tStarts)-1:
        try:
            exonStart, exonEnd = tStarts[h], tStarts[h]+blockSizes[h]
            nextExonStart, nextExonEnd = tStarts[h+1], tStarts[h+1]+blockSizes[h+1]
            jStart, jEnd = exonEnd, nextExonStart
            junctions.add((tName, jStart, jEnd))
        except IndexError:
            print tStarts, blockSizes
            break
        h += 1

def parse(fileName, junctions):
    with open(fileName) as f:
        for alnObj in psl_parser.read(f, 'track'):
            tStarts = alnObj.attrib['tStarts']
            blockSizes = alnObj.attrib['blockSizes']
            tName = alnObj.attrib['tName']
            construct(tName, tStarts, blockSizes, junctions)

firstSample = set([])
secondSample = set([])

parse(sys.argv[1], firstSample)
parse(sys.argv[2], secondSample)

with open('Intersection.out', 'w') as fp:
    for j in firstSample.intersection(secondSample):
        print >> fp, '%s\t%d\t%d' % j

with open('OneDiffTwo.out', 'w') as fp:
    for j in firstSample.difference(secondSample):
        print >> fp, '%s\t%d\t%d' % j

with open('TwoDiffOne.out', 'w') as fp:
    for j in secondSample.difference(firstSample):
        print >> fp, '%s\t%d\t%d' % j
