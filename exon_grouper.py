#! /usr/local/bin/python

'''
    This script groups all exons from the same transcript together.
    The output is in BED format which can be visualized in UCSC genome browser.
    The script requires the alignment of transcript assembly from
    velvet + oases to the referecne genome.
    The alignment has to be in PSL format from GMAP, BLAT. 

    Usage: python exon_grouper.py [transcripts.psl]
    The output is written in a standard output.

    The script is written in Python 2.6.6

    Author: Likit Preeyanon
    Email: preeyano@msu.edu
'''

import psl_parser
import sys
import csv
from pygr import seqdb, sequtil

class Exon(object):

    def __init__(self, reference, start, end, junctions, leftTerminal):
        self.reference = reference 
        self.start = start
        self.end = end
        self.junctions = junctions
        self.cluster = []
        self.connectedExons = [(start, end)]
        self.leftTerminal = leftTerminal

def construct(tName, tStarts, blockSizes, exons):
    '''
        Constructs a dictionary containing all exon objects.
    '''

    for i in range(len(tStarts)-1):
        end = tStarts[i] + blockSizes[i]
        start = tStarts[i]
        juncStart = tStarts[i+1]
        juncEnd = tStarts[i+1] + blockSizes[i+1]

        if end in exons and tName == exons[end].reference:

            if start < exons[end].start:
                if exons[end].leftTerminal:
                    exons[end].start = start

            if start > exons[end].start:
                if exons[end].leftTerminal and i > 0:
                    exons[end].leftTerminal = False
                    exons[end].start = start
                elif not exons[end].leftTerminal and i > 0:
                    exons[end].start = start

            if (juncStart, juncEnd) not in exons[end].junctions:
                exons[end].junctions.append((juncStart, juncEnd))
        else:
            if i == 0:
                leftTerminal = True
            else:
                leftTerminal = False
            junc = [(juncStart, juncEnd)]
            exons[end] = Exon(tName, start, end, junc, leftTerminal)

    lastExonStart = tStarts[-1]
    lastExonEnd = lastExonStart + blockSizes[-1]

    if lastExonEnd in exons and tName == exons[lastExonEnd].reference:
        if lastExonStart > exons[lastExonEnd].start:
            exons[lastExonEnd].start = lastExonStart

    else:
        exons[lastExonEnd] = Exon(tName, \
                            lastExonStart, lastExonEnd, [], False)

def join(exons, exonEnd, groupedExons, newCluster, \
            allConnectedExons, addClusters):
    '''
        This function walks through all exons connected to the starting
        exon.
        Parameters:
            exons: dictionary containing Exon object.
            exonEnd: the end of the starting exon.
            groupedExons: a list containing exons that already clustered.
            newCluster: the end of a starting exon to be used as a
            cluster name.
            allConnectedExons: a list to be added all exons connected to
            the starting exon.
            addClusters: a list to be added all clusters that each exon
            belongs to. 
    '''

    if exonEnd not in groupedExons:
        groupedExons.append(exonEnd)

    if exons[exonEnd].junctions:
        for start, end in exons[exonEnd].junctions:
            allConnectedExons.append((start, end))
            if newCluster not in exons[end].cluster:
                exons[end].cluster.append(newCluster)  # add new cluster to the exon's cluster list.
            for c in exons[end].cluster:
                if c not in addClusters:
                    addClusters.append(c)
            join(exons, end, groupedExons, newCluster, allConnectedExons, addClusters)
    else:
        return

def cluster(exons):
    '''
        Clusters exons from the same genes, isoforms together.
    '''

    allGroupedExons = {}
    exonClusters = []
    clusterReferences = set()
    n = 0
    for num, e in enumerate(sorted(exons)):
        reference = exons[e].reference
        exonEnd = exons[e].end

        try:
            groupedExons = allGroupedExons[reference]
        except KeyError:
            allGroupedExons[reference] = []
            groupedExons = []

        if e in groupedExons:
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
            join(exons, exons[e].end, groupedExons, e, allConnectedExons, newClusters)

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


    for e in exons:
        if len(exons[e].cluster) == 1 and len(exons[e].connectedExons) >= 2:
            exonClusters.append(e)
            clusterReferences.add(exons[e].reference)

    print >> sys.stderr, 'total clusters =', len(exonClusters)
    return exonClusters, clusterReferences

def buildGeneModels(exons, exonClusters, clusterReferences):

    geneModels = {}
    excluded = []
    for ref in clusterReferences:
        if ref not in geneModels:
            geneModels[ref] = []
        for e in exonClusters:
            if exons[e].reference == ref:
                connectedExons = sorted(exons[e].connectedExons)

                '''
                    Eliminate all junctions that do not exist.
                '''
                '''
                filteredConnectedExons = [juncs for juncs in
                    connectedExons if exons[juncs[-1]].start == juncs[0]]
                '''

                '''
                    Resolve intron retention.
                '''
                '''
                newConnectedExons = []
                k = 0
                while k < len(filteredConnectedExons):
                    ex = filteredConnectedExons[k]
                    if ex in excluded:
                        k += 1
                        continue

                    overlapped = []
                    for j in range(len(filteredConnectedExons)):
                        nextEx = filteredConnectedExons[j]
                        if nextEx == ex or nextEx in excluded:
                            continue
                        if ex[-1] >= nextEx[-1] and \
                            nextEx[0] >= ex[0]: # if ex longer than nextEx
                            overlapped.append(nextEx)
                            excluded.append(nextEx)

                    if len(overlapped) > 1:
                        x = 0
                        while x < len(overlapped):
                            ExStart, ExEnd = overlapped[x]
                            if exons[ExEnd].junctions:
                                newConnectedExons.append((ExStart,
                                ExEnd))
                            else:
                                try:
                                    nextExStart, nextExEnd = overlapped[x+1]
                                    if nextExStart == ExStart:
                                        pass
                                    else:
                                        newConnectedExons.append((ExStart,
                                        ExEnd))
                                except IndexError:
                                    ExEnd = ex[-1]
                                    newConnectedExons.append((ExStart,
                                    ExEnd))
                            x += 1
                    else:
                        newConnectedExons.append(ex) # no overlapped found
                    k += 1
                '''
                '''
                    Resolve alternative splice sites.
                '''

                newConnectedExons = connectedExons
                cleanedConExons = []
                h = 1
                exStart, exEnd = newConnectedExons[0]
                while h < len(newConnectedExons):
                    nextStart, nextEnd = newConnectedExons[h]
                    if exStart == nextStart:
                        if exons[nextEnd].junctions:
                            exStart, exEnd = nextStart, nextEnd
                        else:
                            if not exons[exEnd].junctions:
                                exStart, exEnd = nextStart, nextEnd
                    else:
                        if exEnd < nextStart:
                            cleanedConExons.append((exStart, exEnd))
                            exStart, exEnd = nextStart, nextEnd 
                        else:
                            if exEnd < nextEnd:
                                exEnd = nextEnd
                    h += 1
                cleanedConExons.append((exStart, exEnd))
                geneModels[ref].append(cleanedConExons)

    return geneModels

def getSequenceExonWise(geneModels, genome):
    for ref in geneModels:
        transcriptNumber = 1
        op = open(ref+'.fasta', 'w')
        for gene in geneModels[ref]:
            exonNumber = 1
            for exon in gene:
                start, end = exon
                seq = genome[ref][start:end]
                exonID = 'gene_%d:exon_%d' % (transcriptNumber,exonNumber)
                sequtil.write_fasta(op,str(seq), id=exonID)
                exonNumber += 1
            transcriptNumber += 1
        op.close()

def printBed(geneModels):

    writer = csv.writer(sys.stdout, dialect='excel-tab')

    for ref in geneModels:
        transcriptNumber = 0
        for m in geneModels[ref]:
            model = {}
            for v in m:
                model[v[0]] = v[-1]
            transcriptNumber += 1
            chromStart = sorted(model)[0]
            blockStarts = [j - chromStart for j in sorted(model)]
            blockSizes = [model[j] - j for j in sorted(model)]

            chromEnd = blockStarts[-1] + blockSizes[-1] + chromStart
            blockCount = len(blockStarts)
            newBlockStarts = [str(i) for i in blockStarts]
            newBlockSizes = [str(i) for i in blockSizes]
            chrom = ref
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

    return newBlockStarts, newBlockSizes

if __name__ == '__main__':

    print >> sys.stderr, 'Constructing exons ...'

    exons = {}

    for alnObj in psl_parser.read(open(sys.argv[1]), 'track'):
        tStarts = alnObj.attrib['tStarts']
        blockSizes = alnObj.attrib['blockSizes']
        tName = alnObj.attrib['tName']
        construct(tName, tStarts, blockSizes, exons)

    print >> sys.stderr, 'total exons = %d' % len(exons)

    print >> sys.stderr, 'Clustering exons ...'
    exonClusters, clusterReferences = cluster(exons)
    print >> sys.stderr, 'Building gene models ...'
    geneModels = buildGeneModels(exons, exonClusters, clusterReferences)
    genome = seqdb.SequenceFileDB(sys.argv[2])
    getSequenceExonWise(geneModels, genome)
    #sizes, starts = printBed(geneModels)
