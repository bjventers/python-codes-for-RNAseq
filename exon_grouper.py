#! /usr/local/bin/python

'''
    This script groups all exons from the same transcript together.
    The output is in BED format which can be visualized in UCSC genome browser.
    The script requires the alignment of transcript assembly from
    velvet + oases output to the referecne genome.
    The input has to be in PSL format from GMAP, BLAT or this script.

    Usage: python exon_grouper.py [transcripts.psl]
    The output is written in a standard output.

    Author: Likit Preeyanon
    Email: preeyano@msu.edu
'''

import psl_parser
import sys
import time
import csv

class Exon(object):

    def __init__(self, ref, start, end, junctions):
        self.ref = ref 
        self.start = start
        self.end = end
        self.junctions = junctions
        self.cluster = []
        self.connectedExons = [(start, end)]

def construct(tName, tStarts, blockSizes, exons):
    '''
        Constructs a dictionary containing all exon objects.
    '''

    for i in range(len(tStarts[:-1])):
        end = tStarts[i] + blockSizes[i]
        start = tStarts[i]

        if end in exons and tName == exons[end].ref:

            if start < exons[end].start: exons[end].start = start

            juncStart = tStarts[i+1]
            juncEnd = tStarts[i+1] + blockSizes[i+1]

            if (juncStart, juncEnd) not in exons[end].junctions:
                exons[end].junctions.append((juncStart, juncEnd))

        else:
            juncStart = tStarts[i+1]
            juncEnd = tStarts[i+1] + blockSizes[i+1]
            junc = [(juncStart, juncEnd)]
            exons[end] = Exon(tName, start, end, junc)

    lastExonStart = tStarts[-1]
    lastExonEnd = lastExonStart + blockSizes[-1]

    if lastExonEnd in exons and tName == \
        exons[lastExonEnd].ref:
        if lastExonStart < exons[lastExonEnd].start:
            exons[lastExonEnd].start = lastExonStart
    else:
        exons[lastExonEnd] = Exon(tName, \
                            lastExonStart, lastExonEnd, [])

def join(exons, exonEnd, grouped, newCluster, \
            allConnectedExons, addClusters):
    '''
        This function walks through all exons connected to the starting
        exon.
        Parameters:
            exons: dictionary containing Exon object.
            exonEnd: key(end) of the starting exon.
            grouped: a list containing exons that already clustered.
            newCluster: a key(end) of a starting exon to be used as a
            cluster name.
            allConnectedExons: a list to be added all exons connected to
            the starting exon.
            addClusters: a list to be added all clusters that each exon
            belongs to. 
    '''

    if exonEnd not in grouped:
        grouped.append(exonEnd)

    if exons[exonEnd].junctions:
        for start, end in exons[exonEnd].junctions:
            allConnectedExons.append((exons[end].start, exons[end].end))
            if newCluster not in exons[end].cluster:
                exons[end].cluster.append(newCluster)  # add new cluster to the exon's cluster list.
            for c in exons[end].cluster:
                if c not in addClusters:
                    addClusters.append(c)
            join(exons, end, grouped, newCluster, allConnectedExons, addClusters)
    else:
        return

def cluster(exons):

    grouped = []
    exonClusters = []
    n = 0
    now = time.time()
    print >> sys.stderr, 'Clustering ...'
    for num, e in enumerate(sorted(exons)):
        if e in grouped:
            if num % 1000 == 0:
                print >> sys.stderr, '...', num,
                print >> sys.stderr, time.time() - now, len(grouped)
                now = time.time()
            continue
        else:
            '''
                All exons that are connected to this exon are stored in
                allConnectedExons.
            '''
            allConnectedExons = []  

            '''
                All connected exons will be put into the cluster that
                this exon belongs to.
                Note that each exon belongs to its own cluster.
                All clusters that each exon belongs to can be looked up in
                cluster attribute of Exon object. 
            '''
            exons[e].cluster.append(e)
            newClusters = []

            '''
                join() walks through all exons that connected together
                and add them in allConnectedExons list.
                Also, all clusters that found in cluster attribute of each
                connected exon will be added to newClusters list.
            '''
            join(exons, exons[e].end, grouped, e, allConnectedExons, newClusters)
            #print >> sys.stderr, exons[e].start, exons[e].end, allConnectedExons

            if newClusters:
                for start, end in allConnectedExons:
                    for nc in newClusters:
                        if nc not in exons[end].cluster:
                            exons[end].cluster.append(nc)
                        if (start, end) not in exons[nc].connectedExons:
                            exons[nc].connectedExons.append((start, end))
                for nc in newClusters:
                    if nc not in exons[e].cluster:
                        exons[e].cluster.append(nc)

                    if (exons[e].start, exons[e].end) not in exons[nc].connectedExons:
                        exons[nc].connectedExons.append((exons[e].start, exons[e].end))

        if num % 1000 == 0:
            print >> sys.stderr, '...', num, time.time() - now, len(grouped)
            now = time.time()

    for e in sorted(exons):
        if len(exons[e].cluster) == 1 and len(exons[e].connectedExons) >= 2:
            exonClusters.append(e)
    print >> sys.stderr, 'total clusters =', len(exonClusters)
    assert len(exonClusters) == 4
    return exonClusters

def printBed(exons, exonClusters):

    transcriptNumber = 0

    writer = csv.writer(sys.stdout, dialect='excel-tab')

    for e in exonClusters:
        newJunctions = {}
        connectedExons = sorted(exons[e].connectedExons)
        if connectedExons:
            for start, end in connectedExons:
                if start not in newJunctions:
                    newJunctions[start] = end
                else:
                    if end > newJunctions[start]:
                        newJunctions[start] = end

        if newJunctions:
            transcriptNumber += 1
            chromStart = connectedExons[0][0]
            blockStarts = [j - chromStart for j in sorted(newJunctions)]
            blockSizes = [newJunctions[j] - j for j in sorted(newJunctions)]
            #blockEnds = [int(new_junctions[j]) for j in sorted(new_junctions)]

            chromEnd = blockStarts[-1] + blockSizes[-1] + chromStart
            blockCount = len(blockStarts)
            newBlockStarts = [str(i) for i in blockStarts]
            newBlockSizes = [str(i) for i in blockSizes]
            chrom = exons[e].ref
            name="%s_%d" % (chrom, transcriptNumber)
            strand = "+"
            score=1000
            itemRgb="0,0,0"
            writer.writerow((chrom,
                            chromStart,
                            chromEnd, 
                            name,
                            score,
                            strand,
                            chromStart, 
                            chromEnd,
                            itemRgb,
                            blockCount,
                            ','.join(newBlockSizes),
                            ','.join(newBlockStarts)))

if __name__ == '__main__':

    print >> sys.stderr, 'Constructing exons ...'

    exons = {}

    for alnObj in psl_parser.read(open(sys.argv[1]), 'track'):
        tStarts = alnObj.attrib['tStarts']
        blockSizes = alnObj.attrib['blockSizes']
        tName = alnObj.attrib['tName']
        construct(tName, tStarts, blockSizes, exons)

    print >> sys.stderr, 'total exons = %d' % len(exons)

    multipleJunctions = [exon for exon in exons.itervalues() if len(exon.junctions) > 1]

    print >> sys.stderr, 'total multiple junctions', len(multipleJunctions)

    exonClusters = cluster(exons)
    printBed(exons, exonClusters)
