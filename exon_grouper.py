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
from pygr import seqdb, sequtil
from operator import itemgetter
from string import maketrans
from Bio.Blast import NCBIWWW, NCBIXML

def deleteGap(tName, tStarts, blockSizes):
    '''
        Delete all small gaps and overlapped exons.
    '''

    exonSet = []
    i = 0
    ref, start ,end = tName, tStarts[0], tStarts[0]+blockSizes[0]
    while i < range(len(tStarts)): 
        try:
            ref, nextStart, nextEnd = tName, tStarts[i+1], tStarts[i+1]+blockSizes[i+1]
        except IndexError:
            exonSet.append((tName, start, end))
            break
        else:
            if nextStart - end < 21:
                end = nextEnd
            else:
                exonSet.append((tName, start, end))
                start, end = nextStart, nextEnd
        i += 1
    return exonSet

def findLongestEnd(allExons, linkedExons, endExons, exonPositions, ignored):
    allExons = sorted(allExons, reverse=True)
    curRef, curStart, curEnd = allExons[0]
    curPos = exonPositions[(curRef, curStart, curEnd)]
    change = []
    i = 1
    while True:
        i += 1
        try:
            nextRef, nextStart, nextEnd = allExons[i]
            nextPos = exonPositions[(nextRef, nextStart, nextEnd)]
        except IndexError:
            break
        else:
            if curStart == nextStart:
                if curPos == 1 and nextPos == 1:
                    change.append((nextRef, nextStart, nextEnd))
            else:
                for c in change:
                    secondLastExon = endExons[c] 
                    linkedExons[secondLastExon].add((curRef, curStart, curEnd))
                    #print c, linkedExons[c], '-->', curRef, curStart, curEnd, '-->', secondLastExon, linkedExons[secondLastExon]

def construct(tName, tStarts, blockSizes, exons,
                clusters, newClusterID, clusterConnections,
                linkedExons, exonPositions, endExons):
    '''
        Constructs a dictionary containing all exon.
    '''
    exonGroup = set([])
    connection = set([])
    exonSet = deleteGap(tName, tStarts, blockSizes)
    if len(exonSet) == 1:
        return newClusterID

    for i in range(len(exonSet)):
        tName, start, end = exonSet[i]
        try:
            tName, juncExonStart, juncExonEnd = exonSet[i+1]
            try:
                linkedExons[(tName, start, end)].add((tName, juncExonStart, juncExonEnd))
            except KeyError:
                linkedExons[(tName, start, end)] = set([(tName, juncExonStart, juncExonEnd)])
        except IndexError:
            try:
                allLinks = linkedExons[(tName, start, end)]
            except KeyError:
                linkedExons[(tName, start, end)] = set([])
                endExons[(tName, start, end)] = exonSet[i-1]

        exonGroup.add((tName, start, end))

        '''
            Assign position of the exons:
                -1 = first exon
                0  = middle exon
                1  = last exon
        '''
        try:
            position = exonPositions[(tName, start, end)]
        except KeyError:
            if i == 0:
                exonPositions[(tName, start, end)] = -1
            elif i == len(exonSet)-1:
                exonPositions[(tName, start, end)] = 1
            else:
                exonPositions[(tName, start, end)] = 0
        else:
            pass
    ''' 
        If at least one exon connects to an existing cluster,
        add all new exons to that cluster and exons database.
    '''

    for tName, start, end in sorted(exonGroup):
        try:
            clusterID = exons[(tName, start)]
        except KeyError:
            pass
        else:
            connection.add(clusterID)

    if connection:
        clusters[clusterID] = clusters[clusterID].union(exonGroup)
        for tName, start, end in sorted(exonGroup):
            exons[(tName, start)] = clusterID
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
                exons[(tName, start)] = (tName, newClusterID)
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

def walkFork(nodes, linkedExons, passed, paths, visited, ignored, txExons):
    nodes = sorted(nodes)
    if nodes:
        while nodes:
            direction = nodes.pop()
            visited.add(direction)
            if direction not in ignored: 
                passed.append(direction)
                walkFork(linkedExons[direction],
                                    linkedExons,
                                    passed, paths,
                                    visited,
                                    ignored,
                                    txExons,
                                    )
                passed.pop()
    else:
        paths.append(passed[:])

def buildPaths(linkedExons, txExons, allPaths, ignored, visited):
    paths = []
    for c in sorted(txExons):
        if c not in visited and c not in ignored:
            passed = [c]
            visited.add(c)
            walkFork(linkedExons[c], linkedExons,
                                        passed, paths,
                                        visited,
                                        ignored,
                                        txExons,
                                        )
    return paths

def buildGeneModels(mergedClusters, exonPositions):

    geneModels = {}
    for cluster in mergedClusters:
        connectedExons = sorted(mergedClusters[cluster], key=itemgetter(1,2))
        cleanedConExons = []
        ref, exonStart, exonEnd = connectedExons[0]
        exonPosition = exonPositions[(ref, exonStart, exonEnd)]
        h = 1
        while h < len(connectedExons):
            ref, nextExonStart, nextExonEnd = connectedExons[h]
            nextExonPosition = exonPositions[(ref, nextExonStart, nextExonEnd)]
            if exonStart == nextExonStart:
                if nextExonPosition > -1:
                    exonStart, exonEnd = nextExonStart, nextExonEnd
                else:
                    if not exonPosition > -1: # exonPosition == -1?
                        exonStart, exonEnd = nextExonStart, nextExonEnd
            else:
                if nextExonStart-exonEnd >=30:
                    cleanedConExons.append((ref, exonStart, exonEnd))
                    exonStart, exonEnd = nextExonStart, nextExonEnd 
                    exonPosition = nextExonPosition
                else:
                    if exonEnd < nextExonEnd:
                        exonEnd = nextExonEnd
                        exonPosition = nextExonPosition

            h += 1
        cleanedConExons.append((ref, exonStart, exonEnd))
        geneModels[cluster] = cleanedConExons

    return geneModels

def getSequenceExonWiseIsoform(allPaths, genome):
    allSequences = {}
    sequences = []
    for cl in allPaths:
        for gene in allPaths[cl]:
            seq = ''
            for exon in gene:
                r, start, end = exon 
                seq += str(genome[r][start:end])
            sequences.append(seq)
            print '>%s:%d:%d\n%s' % (r, start, end, seq)
        allSequences[cl] = sequences
    return allSequences

def getSequenceExonWiseUnigene(allPaths, genome):
    for cl in allPaths:
        ref, ID = cl
        seq = ''
        allExons = sorted(allPaths[cl])
        firstStart = allExons[0][1]
        lastEnd = allExons[-1][-1]
        print '>%s:%d:%d' % (ref, firstStart, lastEnd)
        for exon in allPaths[cl]:
            r, start, end = exon 
            print str(genome[r][start:end])

def printBedUnigene(clusters):

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

def printBedIsoforms(clusters, genome):

    writer = csv.writer(sys.stdout, dialect='excel-tab')
    geneNo = 0
    for ref, ID in clusters:
        geneNo += 1
        if geneNo%1000 == 0:
            print >> sys.stderr, '...', geneNo, 'finished'
        isoformNum = 0
        for isoforms in clusters[(ref,ID)]:
            cl = sorted(isoforms)
            chromStart = cl[0][1]
            blockStarts = [j[1] - chromStart for j in cl]
            blockSizes = [j[2] - j[1] for j in cl]
            chromEnd = blockStarts[-1] + blockSizes[-1] + chromStart

            try:
                frame, startCodon, stopCodon, length = getStartStopCodon(cl, genome)
                if length <= 30:
                    continue
                #print >> sys.stderr, 'startCodon %d, stopCodon %d' % (startCodon, stopCodon)
                #print >> sys.stderr, 'blockSizess', blockSizes, sum(blockSizes)
            except TypeError:
                thickStart = chromStart
                thickEnd = chromEnd
                strand = '.'
                frame = 'NA'
            else:
                if startCodon < blockSizes[0]:
                    thickStart = chromStart + startCodon
                else:
                    for i in range(len(blockSizes)+1):
                        #print >> sys.stderr, startCodon, sum(blockSizes[:i]), i
                        if startCodon > sum(blockSizes[:i]):
                            continue
                        else:
                            startCodon = startCodon - sum(blockSizes[:i-1])
                            thickStart = blockStarts[i-1] + startCodon + chromStart
                            break

                if stopCodon > sum(blockSizes[:-1]):
                    stopCodon = stopCodon - sum(blockSizes[:-1])
                    thickEnd = blockStarts[-1] + stopCodon + chromStart
                else:
                    for i in range(len(blockSizes)):
                        if stopCodon > sum(blockSizes[:i]):
                            continue
                        else:
                            stopCodon = stopCodon - sum(blockSizes[:i-1])
                            thickEnd = blockStarts[i-1] + stopCodon + chromStart
                            break

                strand = "+" if frame > 0 else "-"

            blockCount = len(blockStarts)
            newBlockStarts = [str(i) for i in blockStarts]
            newBlockSizes = [str(i) for i in blockSizes]

            chrom = ref
            name="%s_%d_%d_frame_%s" % (chrom, ID, isoformNum, str(frame))
            score=1000
            itemRgb="0,0,0"
            writer.writerow((chrom,
                            chromStart,
                            chromEnd, 
                            name,
                            score,
                            strand,
                            thickStart, 
                            thickEnd,
                            itemRgb,
                            blockCount,
                            ','.join(newBlockSizes),
                            ','.join(newBlockStarts)))
            isoformNum += 1

def cleanUpLinkedExons(allExons, linkedExons, exonPositions, ignored):
    h = 0
    keys = sorted(allExons, key=itemgetter(2,1)) # sort by End then by Start.
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
                                if nextStart - curStart < 20:
                                    linkedExons[(nextRef, nextStart, nextEnd)] = \
                                    linkedExons[(nextRef, nextStart, nextEnd)].union(curLinkedExons)
                                    ignored.add((curRef, curStart, curEnd))
                                    curRef, curStart, curEnd = keys[h]
                            elif curPosition == 0 and nextPosition == 0:
                                pass
                            else: 
                                linkedExons[(curRef, curStart, curEnd)] = \
                                linkedExons[(curRef, curStart, curEnd)].union(nextLinkedExons)
                                ignored.add((nextRef, nextStart, nextEnd))
                    else:
                        curRef, curStart, curEnd = keys[h]
    h = 0
    keys = sorted(allExons, key=itemgetter(1,2)) # sort by Start then End.
    curRef, curStart, curEnd = keys[h]
    secondLastExons = set([])
    while True:
        h += 1
        #if (curRef, curStart, curEnd) not in ignored:
        try:
            nextRef, nextStart, nextEnd = keys[h]
            nextPosition = exonPositions[(keys[h])]
            nextLinkedExons = linkedExons[(keys[h])]
        except IndexError:
            break
        else:
            if curRef == nextRef:
                if curStart == nextStart:
                    if curEnd < nextEnd:
                        curPosition = exonPositions[(curRef, curStart, curEnd)]
                        curLinkedExons = linkedExons[(curRef, curStart, curEnd)]
                        if nextPosition == 0 and curPosition == 1:
                            linkedExons[(nextRef, nextStart, nextEnd)] = \
                            linkedExons[(nextRef, nextStart, nextEnd)].union(curLinkedExons)
                            ignored.add((curRef, curStart, curEnd))
                            curRef, curStart, curEnd = keys[h]
                        elif nextPosition == 1 and curPosition == 0:
                            if nextEnd - curEnd < 20:
                                linkedExons[(curRef, curStart, curEnd)] = \
                                linkedExons[(curRef, curStart, curEnd)].union(nextLinkedExons)
                                ignored.add((nextRef, nextStart, nextEnd))
                            else:
                                curRef, curStart, curEnd = keys[h]
                        elif nextPosition == 1 and curPosition == 1:
                            ignored.add((curRef, curStart, curEnd))
                            curRef, curStart, curEnd = keys[h]
                else:
                    curRef, curStart, curEnd = keys[h]
    '''
    for k in sorted(keys, key=lambda x: x[-1]):
        print >> sys.stderr, k, linkedExons[k], exonPositions[k],
        if k in ignored:
            print >> sys.stderr, '*'
        else:
            print >> sys.stderr, '\n'
    '''

def getReadingFrameBLAST(seq):

    result = NCBIWWW.qblast('blastx', 'nr', seq)
    blastRecord = NCBIXML.read(result)
    for alignment in blastRecord.alignments:
        for hsp in alignment.hsps:
            print >> sys.stderr, '***alignment****'
            print >> sys.stderr, 'sequence:', alignment.title
            print >> sys.stderr, 'length:', alignment.length
            print >> sys.stderr, 'e value:', hsp.expect
            print >> sys.stderr, hsp.query[0:75] + '...'
            print >> sys.stderr, hsp.match[0:75] + '...'
            print >> sys.stderr, hsp.sbjct[0:75] + '...'
            print >> sys.stderr, hsp.frame
            print >> sys.stderr, hsp.query_start
            print >> sys.stderr, hsp.query_end
        break
    return hsp.frame[0], hsp.query_start, hsp.query_end

def getStartStopCodon(gene, genome):
    seq = ''
    for exon in gene:
        r, start, end = exon 
        seq += str(genome[r][start:end])

    #print >> sys.stderr, 'Doing BLAST (blastx) search against NR...'
    #readingFrame, alignStart, alignEnd = getReadingFrameBLAST(seq)

    seqLengths = []

    '''
        Forward direction.
    '''
    for frame in [0,1,2]:
        i = frame
        start = False
        while i < len(seq):
            codon = seq[i:i+3]
            if not start:
                if codon in ['ATG', 'CTG', 'GTG', 'TTG', 'ATT']:
                    start = True
                    startPos = i
            else:
                if codon in ['TAG', 'TAA', 'TGA']:
                    seqLengths.append((frame+1, startPos, i+3, i-startPos))
                    i = startPos
                    start = False
            i += 3

    '''
        Reverse direction.
    
    '''
    complement = maketrans('ACGT', 'TGCA')
    revSeq = list(seq)
    revSeq.reverse()
    revSeq = ''.join(revSeq)
    revSeq = revSeq.translate(complement)

    for frame in [0, 1, 2]:
        i = frame
        start = False
        while True:
            codon = revSeq[i:i+3]
            if len(codon) < 3:
                break
            if not start:
                if codon in ['ATG', 'CTG', 'GTG', 'TTG', 'ATT']:
                    start = True
                    startPos = i
            else:
                if codon in ['TAG', 'TAA', 'TGA']:
                    seqLengths.append(((frame + 1)*-1, len(seq) - (i+3), len(seq)-startPos, i-startPos))
                    i = startPos
                    start = False
            i += 3

    if seqLengths:
        frame, startCodon, stopCodon, length = sorted(seqLengths, key=lambda x: x[-1])[-1]
        return sorted(seqLengths, key=lambda x: x[-1])[-1]
    else:
        return None

def main():
    exons = {}
    clusters = {}
    newClusterID = 0
    clusterConnections = {}
    linkedExons = {}
    exonPositions = {}
    endExons = {}
    print >> sys.stderr, 'Parsing and clustering exons..'
    for alnObj in psl_parser.read(open(sys.argv[1]), 'track'):
        tStarts = alnObj.attrib['tStarts']
        blockSizes = alnObj.attrib['blockSizes']
        tName = alnObj.attrib['tName']
        qName = alnObj.attrib['qName']
        newClusterID = construct(tName, tStarts, blockSizes,
                                exons, clusters, newClusterID,
                                clusterConnections,
                                linkedExons, exonPositions,
                                endExons)
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
    ignored = set([])
    for cl in mergedClusters:
        allExons = mergedClusters[cl]
        cleanUpLinkedExons(allExons, linkedExons, exonPositions, ignored)

    print >> sys.stderr, 'Modifying the right end of each transcript..'
    for cl in mergedClusters:
        findLongestEnd(mergedClusters[cl], linkedExons, endExons, exonPositions, ignored)
    print >> sys.stderr, '\nConstructing transcripts..'
    allPaths = {}
    gtTenIsoform = 0
    visited = set([])
    for n, cl in enumerate(mergedClusters):
        txExons = sorted(mergedClusters[cl])
        paths = buildPaths(linkedExons, txExons, allPaths, ignored, visited)
        if len(paths) > 10:
            gtTenIsoform += 1
        allPaths[cl] = paths
        if n %1000 == 0:
            if n > 0:
                print >> sys.stderr, '... %d built..' % n

    geneModels = buildGeneModels(mergedClusters, exonPositions)

    print >> sys.stderr, '\nFinding open reading frames..\n'
    #print >> sys.stderr, '\nWriting gene models in BED format..\n'
    #print >> sys.stderr, mergedClusters
    #printBedUnigene(geneModels)
    #getSequenceExonWiseUnigene(geneModels, genome)
    genome = seqdb.SequenceFileDB(sys.argv[2], verbose=False)
    printBedIsoforms(allPaths, genome)
    #allSequences = getSequenceExonWiseIsoform(allPaths, genome)
    #checkReadingFrame(allSequences)
    print >> sys.stderr, '> 10 isoforms = ', gtTenIsoform

if __name__ == '__main__':
    main()
