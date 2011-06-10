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
        end = tStarts[i] + blockSizes[i]
        start = tStarts[i]
        exonGroup.add((tName, start, end))

    ''' 
        If at least one exon connects to an existing cluster,
        add all new exons to that cluster and exons database.
    '''
    for tName, start, end in sorted(exonGroup):
        #print tName, start, end,
        try:
            clusterID = exons[(tName, end)]
        except KeyError:
            pass

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
    #print '='*40

    return newClusterID

def mergeCluster(clusters, clusterConnections):
    #print '=|='*40
    #print 'Merging...'
    #print '=|='*40
    visited = []
    mergedClusters = {}
    for cluster in clusterConnections:
        if cluster not in visited:
            try:
                linkedExons = mergedClusters[cluster]
            except KeyError:
                mergedClusters[cluster] = clusters[cluster]
            #print '-->', cluster
            for linked in clusterConnections[cluster]:
                #print '\t\t-->', linked
                visited.append(linked)
                mergedClusters[cluster].union(clusters[linked])
    return mergedClusters

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
                print '---newConnectedExons---'
                for e in newConnectedExons:
                    print exons[e].start, exons[e].end
                cleanedConExons = []
                h = 1
                exEnd = newConnectedExons[0]
                exStart = exons[exEnd].start
                print exStart, exEnd
                while h < len(newConnectedExons):
                    nextEnd = newConnectedExons[h]
                    nextStart = exons[nextEnd].start
                    print '\t-->', nextStart, nextEnd,
                    if exStart == nextStart:
                        if exons[nextEnd].junctions:
                            exStart, exEnd = nextStart, nextEnd
                        else:
                            if not exons[exEnd].junctions:
                                exStart, exEnd = nextStart, nextEnd
                        print
                    else:
                        if nextStart-exEnd >=30:
                            cleanedConExons.append((exStart, exEnd))
                            print nextStart-exEnd, 'real exon', exStart, exEnd, 'added'
                            exStart, exEnd = nextStart, nextEnd 
                        else:
                            if exEnd < nextEnd:
                                print nextStart-exEnd,
                                exEnd = nextEnd
                                print 'new exEnd', exStart, exEnd
                    h += 1
                    print exStart, exEnd
                cleanedConExons.append((exStart, exEnd))
                geneModels[ref].append(cleanedConExons)

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
def printBed2(clusters):

    writer = csv.writer(sys.stdout, dialect='excel-tab')
    for ref, ID in clusters:
        cl = sorted(clusters[(ref,ID)])
        transcriptNumber = ID
        chromStart = cl[0][1]
        blockStarts = [j[1] - chromStart for j in cl]
        blockSizes = [j[2] - j[1] for j in cl]

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
    print >> sys.stdout, 'Parsing and clustering exons..'
    for alnObj in psl_parser.read(open(sys.argv[1]), 'track'):
        tStarts = alnObj.attrib['tStarts']
        blockSizes = alnObj.attrib['blockSizes']
        tName = alnObj.attrib['tName']
        qName = alnObj.attrib['qName']
        newClusterID = construct(tName, tStarts, blockSizes,
                                    exons, clusters, newClusterID,
                                    clusterConnections,
                                    )
    print >> sys.stdout, '%d exons have been loaded..' % len(exons)
    print >> sys.stdout, 'Linking clusters...'
    mergedClusters = mergeCluster(clusters, clusterConnections)
    '''
    print '--|--'*20
    for mc in mergedClusters:
        print mc, mergedClusters[mc]
        print '-*-'*40
    c0 = clusters[('chr1', 0)]
    c1 = clusters[('chr1', 1)]
    c2 = clusters[('chr1', 2)]
    print 'c1 vs c2'
    print c1.intersection(c2)
    print 'c2 vs c0'
    print c2.intersection(c0)
    print 'c1 vs c0'
    print c0.intersection(c1)
    #print >> sys.stderr, 'total exons = %d' % len(exons)
    #print >> sys.stderr, 'Building gene models ...'
    #geneModels = buildGeneModels(exons, exonClusters, clusterReferences)
    #validateExonLength(geneModels)
    #genome = seqdb.SequenceFileDB(sys.argv[2], verbose=False)
    #getSequenceExonWise2(geneModels, genome)
    '''
    sizes, starts = printBed2(mergedClusters)
    print >> sys.stderr, 'total %d transcripts found.' % (len(mergedClusters))
