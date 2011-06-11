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

def construct(tName, tStarts, blockSizes, exons,
                clusters, newClusterID, clusterConnections):
    '''
        Constructs a dictionary containing all exon objects.
    '''
    exonGroup = set([])
    connection = set([])

    for i in range(len(tStarts)):
        if i == 0:
            exonPosition = 0
        elif i == len(tStarts) - 1:
            exonPosition = -1
        else:
            exonPosition = 1

        end = tStarts[i] + blockSizes[i]
        start = tStarts[i]
        exonGroup.add((tName, start, end, exonPosition))

    ''' 
        If at least one exon connects to an existing cluster,
        add all new exons to that cluster and exons database.
    '''
    for tName, start, end, exonPosition in sorted(exonGroup):
        try:
            clusterID = exons[(tName, end)]
        except KeyError:
            pass

    for tName, start, end, exonPosition in sorted(exonGroup):
        try:
            clusterID = exons[(tName, end)]
        except KeyError:
            pass
        else:
            connection.add(clusterID)

    if connection:
        clusters[clusterID] = clusters[clusterID].union(exonGroup)
        for tName, start, end, exonPosition in sorted(exonGroup):
            exons[(tName, end)] = clusterID
        for c in connection:
            clusterConnections[c] = clusterConnections[c].union(connection)
    else:
        '''
            If no exons connects to any exon in existing clusters,
            new cluster will be created for the new group of exons.
            All new exons are also added to the exons database.
        '''
        clusters[(tName, newClusterID)] = exonGroup
        clusterConnections[(tName, newClusterID)] = set([])
        for tName, start, end, exonPosition in sorted(exonGroup):
            try:
                exons[(tName, end)] = (tName, newClusterID)
            except:
                raise KeyError
        newClusterID += 1

    return newClusterID

def mergeCluster(clusters, clusterConnections):
    visited = []
    mergedClusters = {}
    for cluster in clusterConnections:
        if cluster not in visited:
            try:
                linkedExons = mergedClusters[cluster]
            except KeyError:
                mergedClusters[cluster] = clusters[cluster]
            for linked in clusterConnections[cluster]:
                visited.append(linked)
                mergedClusters[cluster].union(clusters[linked])
    return mergedClusters

def buildGeneModels(mergedClusters):

    geneModels = {}
    for cluster in mergedClusters:
        connectedExons = sorted(mergedClusters[cluster])
        cleanedConExons = []
        ref, exonStart, exonEnd, exonPosition = connectedExons[0]
        h = 1
        while h < len(connectedExons):
            nextRef, nextExonStart, nextExonEnd, nextExonPosition = connectedExons[h]
            if exonStart == nextExonStart:
                if nextExonPosition > -1:
                    exonStart, exonEnd = nextExonStart, nextExonEnd
                else:
                    if not exonPosition > -1:
                        exonStart, exonEnd = nextExonStart, nextExonEnd
            else:
                if nextExonStart-exonEnd >=30:
                    cleanedConExons.append((exonStart, exonEnd))
                    exonStart, exonEnd = nextExonStart, nextExonEnd 
                else:
                    if exonEnd < nextExonEnd:
                        exonEnd = nextExonEnd
            h += 1
        cleanedConExons.append((exonStart, exonEnd))
        geneModels[cluster] = cleanedConExons

    return geneModels

def validateExonLength(geneModels):
    for ref in geneModels:
        geneNo = 1
        for gene in geneModels[ref]:
            geneLength = 0
            for i in range(len(gene)-1):
                if i == 0:
                    continue
                else:
                    start, end = gene[i]
                    geneLength = (end - start)
            print '%s gene %d length = %d' % (ref, geneNo, geneLength%3)
            geneNo += 1

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

def getSequenceExonWise2(geneModels, genome):
    for ref in geneModels:
        transcriptNumber = 1
        op = open(ref+'.fasta', 'w')
        for gene in geneModels[ref]:
            exonSeqs = []
            shortExon = False
            for exon in gene:
                start, end = exon
                seq = genome[ref][start:end]
                exonSeqs.append(str(seq))
                if len(seq)<3:
                    print transcriptNumber, gene
                    '''
                    if not shortExon:
                        shortExon = True
                    '''
            '''
            if shortExon:
                print >> op, '>gene_%d' % transcriptNumber
                for seq in exonSeqs:
                    print >> op, str(seq) 
            '''
            transcriptNumber += 1
        op.close()

def printBed(clusters):

    writer = csv.writer(sys.stdout, dialect='excel-tab')
    for ref, ID in clusters:
        cl = sorted(clusters[(ref,ID)])
        transcriptNumber = ID
        chromStart = cl[0][0]
        blockStarts = [j[0] - chromStart for j in cl]
        blockSizes = [j[1] - j[0] for j in cl]

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

    exons = {}
    clusters = {}
    newClusterID = 0
    clusterConnections = {}
    print >> sys.stderr, 'Parsing and clustering exons..'
    for alnObj in psl_parser.read(open(sys.argv[1]), 'track'):
        tStarts = alnObj.attrib['tStarts']
        blockSizes = alnObj.attrib['blockSizes']
        tName = alnObj.attrib['tName']
        qName = alnObj.attrib['qName']
        newClusterID = construct(tName, tStarts, blockSizes,
                                    exons, clusters, newClusterID,
                                    clusterConnections,
                                    )
    sumExons = {}
    for ref, end in exons:
        try:
            sumExons[ref] += 1
        except KeyError:
            sumExons[ref] = 1
    for ref in sorted(sumExons):
        print >> sys.stderr, '\t%s has %d exon(s).' % (ref, sumExons[ref])

    print >> sys.stderr, '\nLinking clusters..'
    mergedClusters = mergeCluster(clusters, clusterConnections)
    print >> sys.stderr, '\nBuilding gene models..'
    geneModels = buildGeneModels(mergedClusters)
    allReferences = {}
    for ref, ID in clusters.keys():
        try:
            allReferences[ref] += 1
        except KeyError:
            allReferences[ref] = 1
    for ref in sorted(allReferences):
        print >> sys.stderr, '\t%s has %d genes.' % (ref, allReferences[ref])
    #validateExonLength(geneModels)
    #genome = seqdb.SequenceFileDB(sys.argv[2], verbose=False)
    #getSequenceExonWise2(geneModels, genome)
    print >> sys.stderr, '\nWriting gene models in BED format..\n'
    sizes, starts = printBed(geneModels)
