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
    Last update: 13/11/2011
'''

import psl_parser
import sys
import csv
#from pygr import seqdb, sequtil
from operator import itemgetter

def construct(tName, tStarts, blockSizes, exons,
                clusters, newClusterID, clusterConnections,
                linkedExons, exonPositions):
    '''
        Constructs a dictionary containing all exon objects.
    '''
    exonGroup = set([])
    connection = set([])

    for i in range(len(tStarts)):
        end = tStarts[i] + blockSizes[i]
        start = tStarts[i]
        try:
            juncExonEnd = tStarts[i+1]+blockSizes[i+1]
            juncExonStart = tStarts[i+1]
            try:
                linkedExons[(tName, start, end)].add((tName, juncExonStart, juncExonEnd))
            except KeyError:
                linkedExons[(tName, start, end)] = set([])
        except IndexError:
            try:
                allLinks = linkedExons[(tName, start, end)]
            except KeyError:
                linkedExons[(tName, start, end)] = set([])

        exonGroup.add((tName, start, end))
        try:
            position = exonPositions[(tName, start, end)]
        except KeyError:
            if i == 0:
                exonPositions[(tName, start, end)] = -1
            elif i == len(tStarts)-1:
                exonPositions[(tName, start, end)] = 1
            else:
                exonPositions[(tName, start, end)] = 0

    ''' 
        If at least one exon connects to an existing cluster,
        add all new exons to that cluster and exons database.
    '''

    for tName, start, end in sorted(exonGroup):
        try:
            clusterID = exons[(tName, end)]
        except KeyError:
            pass
        else:
            connection.add(clusterID)

    if connection:
        clusters[clusterID] = clusters[clusterID].union(exonGroup)
        for tName, start, end in sorted(exonGroup):
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
        for tName, start, end in sorted(exonGroup):
            try:
                exons[(tName, end)] = (tName, newClusterID)
            except:
                raise KeyError
        newClusterID += 1

    return newClusterID

def walk(allExons, nodes, clusters, clusterConnections, visited):

    ''' This function walks over all clusters connected to
    the starting cluster and combine all exons to those is
    the starting cluster.

    '''
    if nodes not in visited:
        for n in nodes:
            if n not in visited:
                visited.append(n)
                allExons = allExons.union(clusters[n])
                nodes = clusterConnections[n]
                allExons = walk(allExons, nodes, clusters, clusterConnections, visited)
    return allExons

def mergeClusters(clusters, clusterConnections):

    ''' This function merge all clusters sharing at least one exon
    together to form a bigger cluster.

    '''
    visited = []
    mergedClusters = {}
    for i, c in enumerate(clusterConnections, start=1):
        if c not in visited:
            visited.append(c)
            allExons = clusters[c]
            nodes = clusterConnections[c]
            allExons = walk(allExons, nodes, clusters, clusterConnections, visited) 
            mergedClusters[c] = allExons
        if i%1000 == 0:
            print >> sys.stderr, '...', i, 'merged..'
    return mergedClusters

def walkFork(nodes, linkedExons, passed, paths, visited):
    nodes = sorted(nodes)
    if nodes:
        while nodes:
            direction = nodes.pop()
            visited.append(direction)
            passed.append(direction)
            walkFork(linkedExons[direction],
                                linkedExons,
                                passed, paths,
                                visited,
                                )
            passed.pop()
    else:
        paths.append(passed[:])

def buildPaths(linkedExons, txExons, allPaths):
    visited = []
    paths = []
    for c in sorted(txExons):
        if c not in visited and c in linkedExons:
            passed = [c]
            visited.append(c)
            walkFork(linkedExons[c], linkedExons,
                                        passed, paths,
                                        visited,
                                        )
    return paths

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
        isoformNum = 0
        for isoforms in clusters[(ref,ID)]:
            cl = sorted(isoforms)
            chromStart = cl[0][1]
            blockStarts = [j[1] - chromStart for j in cl]
            blockSizes = [j[2] - j[1] for j in cl]

            chromEnd = blockStarts[-1] + blockSizes[-1] + chromStart
            blockCount = len(blockStarts)
            newBlockStarts = [str(i) for i in blockStarts]
            newBlockSizes = [str(i) for i in blockSizes]
            chrom = ref
            name="%s_%d_%d" % (chrom, ID, isoformNum)
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
            isoformNum += 1

    return newBlockStarts, newBlockSizes

def cleanUpLinkedExons(linkedExons, exonPositions):
    ignored = []
    h = 0
    keys = sorted(linkedExons.keys(), key=itemgetter(2,1)) # sort on Start then End.
    curRef, curStart, curEnd = keys[0]
    while True:
        h += 1
        if (curRef, curStart, curEnd) not in ignored:
            try:
                nextRef, nextStart, nextEnd = keys[h]
                nextPosition = exonPositions[(keys[h])]
                nextLinkedExons = linkedExons[(keys[h])]
            except IndexError:
                break
            else:
                if curRef == nextRef:
                    if curEnd == nextEnd:
                        if curStart < nextStart:
                            curPosition = exonPositions[(curRef, curStart, curEnd)]
                            curLinkedExons = linkedExons[(curRef, curStart, curEnd)]
                            if nextPosition == 0 and curPosition == -1:
                                if nextStart - curStart < 10:
                                    linkedExons[(nextRef, nextStart, nextEnd)] = \
                                    linkedExons[(nextRef, nextStart, nextEnd)].union(curLinkedExons)
                                    ignored.append((curRef, curStart, curEnd))
                                    curRef, curStart, curEnd = keys[h]
                            else:
                                linkedExons[(curRef, curStart, curEnd)] = \
                                linkedExons[(curRef, curStart, curEnd)].union(nextLinkedExons)
                                ignored.append((nextRef, nextStart, nextEnd))
                    else:
                        curRef, curStart, curEnd = keys[h]
        print >> sys.stderr, h, 'cleaned..'
    cleanedLinkedExons = {}
    for l in keys:
        if l not in ignored:
            cleanedLinkedExons[l] = linkedExons[l]
    '''
    for l in sorted(cleanedLinkedExons):
        print l, cleanedLinkedExons[l], exonPositions[l]
    '''
    return cleanedLinkedExons

def main():
    exons = {}
    clusters = {}
    newClusterID = 0
    clusterConnections = {}
    linkedExons = {}
    exonPositions = {}
    print >> sys.stderr, 'Parsing and clustering exons..'
    for alnObj in psl_parser.read(open(sys.argv[1]), 'track'):
        tStarts = alnObj.attrib['tStarts']
        blockSizes = alnObj.attrib['blockSizes']
        tName = alnObj.attrib['tName']
        qName = alnObj.attrib['qName']
        newClusterID = construct(tName, tStarts, blockSizes,
                                    exons, clusters, newClusterID,
                                    clusterConnections,
                                    linkedExons, exonPositions)
    sumExons = {}
    for ref, end in exons:
        try:
            sumExons[ref] += 1
        except KeyError:
            sumExons[ref] = 1
    for ref in sorted(sumExons):
        print >> sys.stderr, '\t%s has %d exon(s).' % (ref, sumExons[ref])
    print >> sys.stderr, '\nTotal %d cluster(s) found.' % len(clusters)
    print >> sys.stderr, '\nMerging clusters..'
    mergedClusters = mergeClusters(clusters, clusterConnections)
    print >> sys.stderr, '\nCleaning up..'
    linkedExons = cleanUpLinkedExons(linkedExons, exonPositions)
    print >> sys.stderr, '\nBuilding gene models..'
    allPaths = {}
    for n, cl in enumerate(mergedClusters):
        txExons = mergedClusters[cl]
        paths = buildPaths(linkedExons, txExons, allPaths)
        allPaths[cl] = paths
        print >> sys.stderr, n, 'built..'
    '''
    print >> sys.stderr, '\nBuilding gene models..'
    geneModels = buildGeneModels(mergedClusters)
    allReferences = {}
    for ref, ID in clusters.keys():
        try:
            allReferences[ref] += 1
        except KeyError:
            allReferences[ref] = 1
    for ref in sorted(allReferences):
        print >> sys.stderr, '\t%d gene(s) found in %s.' % (allReferences[ref], ref)
    #validateExonLength(geneModels)
    #genome = seqdb.SequenceFileDB(sys.argv[2], verbose=False)
    print >> sys.stderr, '\nWriting gene models in BED format..\n'
    '''
    getSequenceExonWise(allPaths, genome)
    #sizes, starts = printBed(allPaths)

if __name__ == '__main__':
    main()
