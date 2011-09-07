#! /usr/local/bin/python

''' This script groups all exons from the same transcript together.

The output is in BED format which can be visualized in UCSC genome browser.
The script requires the alignment of transcript assembly from
velvet + oases to the referecne genome.
The alignment has to be in PSL format from GMAP, BLAT.

Usage: python exon_grouper.py [transcripts.psl]
The output is written in a standard output.

The script is tested in Python 2.7.2

Author: Likit Preeyanon
Email: preeyano@msu.edu

'''

import sys
import csv
import sqlite3
from optparse import OptionParser
from operator import itemgetter

from pygr import seqdb
from Bio import SeqIO
from Bio.Blast import NCBIWWW, NCBIXML
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC
from Bio.Data import CodonTable
from Bio.SeqRecord import SeqRecord

import psl_parser

'''Setup option parser'''
parser = OptionParser()
parser.add_option('-g', '--genome', dest='genome',
                    help='genome sequence', metavar='FILE')
parser.add_option('-n', '--basename', dest='basename',
                    help='basename for all output files')
parser.add_option('-i', '--input', dest='infile',
                    help='input file', metavar='FILE')
parser.add_option('-u', '--MINIMUM_UTR_LENGTH', dest='minimumUTRLength',
                    type='int',
                    help='minimum number of bases to be considerd UTR.')


class Isoform(object):
    """Isoform object"""
    def __init__(self, chrom, geneID, isoformID, exons, genome):
        self.chrom = chrom
        self.geneID = geneID
        self.isoformID = isoformID
        self.exons = sorted(exons)
        self.chromStart = exons[0][1]
        self.chromEnd = exons[-1][-1]
        self.frame, self.startCodon, self.stopCodon, \
                self.length, self.dnaSeq = self.__getStartStopCodon(genome)
        self.redundant = False

        if self.frame:
            self.mrnaSeq = self.dnaSeq.transcribe()
            self.orf, self.proteinSeq = self.__getProteinSeq()
        else:
            self.mrnaSeq = None
            self.orf = None
            self.proteinSeq = None

        self.strand = '+' if self.frame > 0 else '-'

    def __getStartStopCodon(self, genome):
        seq = ''
        for exon in self.exons:
            r, start, end = exon
            seq += str(genome[r][start:end])

        '''Create biopython seqObject'''
        bioSeq = Seq(seq.upper(), IUPAC.ambiguous_dna)
        seqLengths = []

        '''Define Standard Start/Stop codons'''
        standardTable = CodonTable.unambiguous_dna_by_name['Standard']

        ''' Forward direction.'''
        for frame in range(3):
            i = frame
            start = False
            while True:
                codon = str(bioSeq[i:i + 3])
                if len(codon) < 3:
                    break
                if not start:
                    if codon in standardTable.start_codons:
                        start = True
                        startPos = i
                else:
                    if codon in standardTable.stop_codons:
                        seqLengths.append((frame + 1, startPos,
                                            i + 3, i - startPos))
                        i = startPos
                        start = False
                i += 3

        ''' Reverse direction'''

        '''Get a reverse complement of bioseq'''
        bioRevSeq = bioSeq.reverse_complement()

        for frame in range(3):
            i = frame
            start = False
            while True:
                codon = str(bioRevSeq[i:i + 3])
                if len(codon) < 3:
                    break
                if not start:
                    if codon in standardTable.start_codons:
                        start = True
                        startPos = i
                else:
                    if codon in standardTable.stop_codons:
                        seqLengths.append(((frame + 1) * -1,
                                            startPos, i + 3,
                                            i - startPos))
                        i = startPos
                        start = False
                i += 3

        if seqLengths:
            frame, startCodon, stopCodon, length \
                    = sorted(seqLengths, key=lambda x: x[-1])[-1]
            if frame > 0:
                return frame, startCodon, stopCodon, length, bioSeq
            else:
                return frame, startCodon, stopCodon, length, bioRevSeq
        else:
            return (None, None, None, None, bioSeq)

    def __getProteinSeq(self):
        if self.frame:
            orf = self.mrnaSeq[self.startCodon:self.stopCodon]
            return orf, orf.translate(cds=True)
        else:
            return None

    def getReferenceBasedStartStopCodon(self):
        start = len(self.mrnaSeq) - self.stopCodon
        end = len(self.mrnaSeq) - self.startCodon
        return start, end


def deleteGap(tName, tStarts, blockSizes):
    '''
        Delete all small gaps and overlapped exons.
    '''

    exonSet = []
    i = 0
    ref, start, end = tName, tStarts[0], tStarts[0] + blockSizes[0]

    while i < range(len(tStarts)):
        try:
            ref, nextStart, nextEnd = tName, tStarts[i + 1], \
                                        tStarts[i + 1] + blockSizes[i + 1]
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
            tName, juncExonStart, juncExonEnd = exonSet[i + 1]
            try:
                linkedExons[(tName, start, end)].add((tName,
                                                    juncExonStart,
                                                    juncExonEnd))
            except KeyError:
                linkedExons[(tName, start, end)] = set([(tName,
                                                        juncExonStart,
                                                        juncExonEnd)])
        except IndexError:
            try:
                allLinks = linkedExons[(tName, start, end)]
            except KeyError:
                linkedExons[(tName, start, end)] = set([])
                endExons[(tName, start, end)] = exonSet[i - 1]

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
            elif i == len(exonSet) - 1:
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
                allExons = walk(allExons, nodes,
                            clusters, clusterConnections,
                            visited)
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
            allExons = walk(allExons, nodes,
                            clusters, clusterConnections,
                            visited)
            mergedClusters[c] = allExons
        if i % 1000 == 0:
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
        connectedExons = sorted(mergedClusters[cluster], key=itemgetter(1, 2))
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
                    if not exonPosition > -1:  # exonPosition == -1?
                        exonStart, exonEnd = nextExonStart, nextExonEnd
            else:
                if nextExonStart - exonEnd >= 30:
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
        cl = sorted(clusters[(ref, ID)])
        transcriptNumber = ID
        chromStart = cl[0][1]
        blockStarts = [j[1] - chromStart for j in cl]
        blockSizes = [j[2] - j[1] for j in cl]

        chromEnd = blockStarts[-1] + blockSizes[-1] + chromStart
        blockCount = len(blockStarts)
        newBlockStarts = [str(i) for i in blockStarts]
        newBlockSizes = [str(i) for i in blockSizes]
        chrom = ref
        name = "%s_%d" % (chrom, transcriptNumber)
        strand = "+"
        score = 1000
        itemRgb = "0,0,0"
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


def writeBEDFile(allGenes, basename):
    writer = csv.writer(open(basename + '.models.bed', 'w'),
                        dialect='excel-tab')
    for chrom in allGenes:
        for geneID in allGenes[chrom]:
            for isoform in allGenes[chrom][geneID]:
                if isoform.redundant:
                    continue

                blockStarts = [j[1] - isoform.chromStart \
                                    for j in isoform.exons]
                blockSizes = [j[2] - j[1] for j in isoform.exons]

                if isoform.frame:
                    if isoform.strand == '+':
                        startCodon, stopCodon = isoform.startCodon, \
                                                isoform.stopCodon
                    else:
                        startCodon, stopCodon = \
                            isoform.getReferenceBasedStartStopCodon()

                    if startCodon < blockSizes[0]:
                        thickStart = isoform.chromStart + startCodon
                    else:
                        for i in range(len(blockSizes) + 1):
                            if startCodon > sum(blockSizes[:i]):
                                continue
                            else:
                                newStartCodon = startCodon - \
                                        sum(blockSizes[:i - 1])
                                thickStart = blockStarts[i - 1] \
                                        + newStartCodon + isoform.chromStart
                                break

                    if stopCodon > sum(blockSizes[:-1]):
                        newStopCodon = stopCodon - sum(blockSizes[:-1])
                        thickEnd = blockStarts[-1] + \
                                    newStopCodon + \
                                    isoform.chromStart
                    else:
                        for i in range(len(blockSizes)):
                            if stopCodon > sum(blockSizes[:i]):
                                continue
                            else:
                                newStopCodon = stopCodon - \
                                                sum(blockSizes[:i - 1])
                                thickEnd = blockStarts[i - 1] + \
                                            newStopCodon + \
                                            isoform.chromStart
                                break
                else:
                    thickStart = isoform.chromStart
                    thickEnd = isoform.chromEnd

                strand = isoform.strand if isoform.frame else '.'
                blockCount = len(blockStarts)
                newBlockStarts = [str(i) for i in blockStarts]
                newBlockSizes = [str(i) for i in blockSizes]

                name = "%s:%d.%d" % (chrom, geneID, isoform.isoformID)
                score = 1000
                itemRgb = "0,0,0"
                writer.writerow((chrom,
                                isoform.chromStart,
                                isoform.chromEnd,
                                name,
                                score,
                                strand,
                                thickStart,
                                thickEnd,
                                itemRgb,
                                blockCount,
                                ','.join(newBlockSizes),
                                ','.join(newBlockStarts)))


def cleanUpLinkedExons(allExons, linkedExons,
                        exonPositions, ignored,
                        minimumUTRLength=100):
    h = 0
    keys = sorted(allExons, key=itemgetter(2, 1))  # sort by End then by Start.
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
                            curPosition = exonPositions[(curRef,
                                                        curStart,
                                                        curEnd)]
                            curLinkedExons = linkedExons[(curRef,
                                                        curStart,
                                                        curEnd)]
                            if nextPosition == 0 and curPosition == -1:
                                if nextStart - curStart < minimumUTRLength:
                                    linkedExons[(nextRef,
                                                nextStart,
                                                nextEnd)] = \
                                    linkedExons[(nextRef,
                                                nextStart,
                                                nextEnd)].union(curLinkedExons)
                                    ignored.add((curRef, curStart, curEnd))
                                    curRef, curStart, curEnd = keys[h]
                            elif curPosition == 0 and nextPosition == 0:
                                pass
                            else:
                                linkedExons[(curRef,
                                            curStart,
                                            curEnd)] = \
                                linkedExons[(curRef,
                                            curStart,
                                            curEnd)].union(nextLinkedExons)
                                ignored.add((nextRef,
                                            nextStart,
                                            nextEnd))
                    else:
                        curRef, curStart, curEnd = keys[h]
    h = 0
    keys = sorted(allExons, key=itemgetter(1, 2))  # sort by Start then End.
    curRef, curStart, curEnd = keys[h]
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
                        curPosition = exonPositions[(curRef,
                                                    curStart,
                                                    curEnd)]
                        curLinkedExons = linkedExons[(curRef,
                                                    curStart,
                                                    curEnd)]
                        if nextPosition == 0 and curPosition == 1:
                            linkedExons[(nextRef,
                                        nextStart,
                                        nextEnd)] = \
                            linkedExons[(nextRef,
                                        nextStart,
                                        nextEnd)].union(curLinkedExons)
                            ignored.add((curRef, curStart, curEnd))
                            curRef, curStart, curEnd = keys[h]
                        elif nextPosition == 1 and curPosition == 0:
                            if nextEnd - curEnd < minimumUTRLength:
                                linkedExons[(curRef,
                                            curStart,
                                            curEnd)] = \
                                linkedExons[(curRef,
                                            curStart,
                                            curEnd)].union(nextLinkedExons)
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
    '''Deprecated'''

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

    '''Create biopython seqObject'''
    bioSeq = Seq(seq, IUPAC.unambiguous_dna)
    seqLengths = []

    '''Define Standard Start/Stop codons'''
    standardTable = CodonTable.unambiguous_dna_by_name['Standard']

    ''' Forward direction.'''
    for frame in range(3):
        i = frame
        start = False
        while True:
            codon = str(bioSeq[i:i + 3])
            if len(codon) < 3:
                break
            if not start:
                if codon in standardTable.start_codons:
                    start = True
                    startPos = i
            else:
                if codon in standardTable.stop_codons:
                    seqLengths.append((frame + 1,
                                        startPos, i + 3,
                                        i - startPos))
                    i = startPos
                    start = False
            i += 3

    ''' Reverse direction'''

    '''Get a reverse complement of bioseq'''
    bioRevSeq = bioSeq.reverse_complement()

    for frame in range(3):
        i = frame
        start = False
        while True:
            codon = str(bioRevSeq[i:i + 3])
            if len(codon) < 3:
                break
            if not start:
                if codon in standardTable.start_codons:
                    start = True
                    startPos = i
            else:
                if codon in standardTable.stop_codons:
                    seqLengths.append(((frame + 1) * -1,
                                        len(seq) - (i + 3),
                                        len(seq) - startPos,
                                        i - startPos))
                    i = startPos
                    start = False
            i += 3

    if seqLengths:
        frame, startCodon, stopCodon, length = \
                sorted(seqLengths, key=lambda x: x[-1])[-1]
        return sorted(seqLengths, key=lambda x: x[-1])[-1]
    else:
        return None


def findRedundantSequence(allGenes):
    for chrom in allGenes:
        for geneID in allGenes[chrom]:
            for isoform1 in allGenes[chrom][geneID]:
                if isoform1.redundant:
                    continue
                else:
                    for isoform2 in allGenes[chrom][geneID]:
                        if isoform1.isoformID == isoform2.isoformID:
                            continue
                        else:
                            if str(isoform1.dnaSeq) == str(isoform2.dnaSeq):
                                isoform2.redundant = True
    return None


def saveToDB(allGenes):
    conn = sqlite3.connect('geneModels')
    cursor = conn.cursor()

    cursor.execute('''
                    create table models (
                    isoformID varchar,
                    chrom varchar,
                    chromStart integer,
                    chromEnd integer,
                    cdsStart integer,
                    cdsEnd integer,
                    strand varchar,
                    exonStarts varchar,
                    exonEnds varchar)
                    ''')
    cursor.execute('''create index isoformID on models (isoformID)''')


def main(options, args):
    exons = {}
    clusters = {}
    newClusterID = 0
    clusterConnections = {}
    linkedExons = {}
    exonPositions = {}
    endExons = {}
    print >> sys.stderr, 'Minimum UTR length = ', options.minimumUTRLength
    print >> sys.stderr, 'Parsing and clustering exons..'
    for n, alnObj in enumerate(psl_parser.read(open(options.infile), 'track')):
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
        cleanUpLinkedExons(allExons,
                            linkedExons,
                            exonPositions,
                            ignored,
                            options.minimumUTRLength)

    print >> sys.stderr, 'Modifying the right end of each transcript..'
    for cl in mergedClusters:
        findLongestEnd(mergedClusters[cl],
                        linkedExons,
                        endExons,
                        exonPositions,
                        ignored)
    print >> sys.stderr, '\nConstructing transcripts..'
    allPaths = {}
    gtTenIsoform = 0
    visited = set([])
    for n, cl in enumerate(mergedClusters):
        txExons = sorted(mergedClusters[cl])
        paths = buildPaths(linkedExons, txExons, allPaths, ignored, visited)
        allPaths[cl] = paths
        if n % 1000 == 0:
            if n > 0:
                print >> sys.stderr, '... %d built..' % n

    #geneModels = buildGeneModels(mergedClusters, exonPositions)
    #print >> sys.stderr, '\nWriting gene models in BED format..\n'
    #print >> sys.stderr, mergedClusters
    #printBedUnigene(geneModels)
    #getSequenceExonWiseUnigene(geneModels, genome)
    genome = seqdb.SequenceFileDB(options.genome, verbose=False)
    #allSequences = getSequenceExonWiseIsoform(allPaths, genome)
    #checkReadingFrame(allSequences)

    '''Create isoform objects from allPaths and
    search for ORF.

    '''
    allGenes = {}
    print >> sys.stderr, 'Searching for ORF...'
    n = 0
    for chrom, geneID in allPaths:
        isoformID = 0
        for isoExons in allPaths[(chrom, geneID)]:
            isoform = Isoform(chrom, geneID, isoformID, isoExons, genome)
            if chrom not in allGenes:
                allGenes[chrom] = {}
                allGenes[chrom][geneID] = [isoform]
            else:
                try:
                    allGenes[chrom][geneID].append(isoform)
                except KeyError:
                    allGenes[chrom][geneID] = [isoform]
            isoformID += 1
        n += 1
        if n % 1000 == 0:
            print >> sys.stderr, '...', n

    print >> sys.stderr, 'Removing redundant sequences..'
    findRedundantSequence(allGenes)

    '''Creating sequence records for each DNA, RNA and protein sequences.'''
    isoformDNASeqs = []
    isoformProteinSeqs = []
    isoformRNASeqs = []
    for chrom in allGenes:
        for geneID in allGenes[chrom]:
            isoformID = 0
            for isoform in allGenes[chrom][geneID]:
                if not isoform.redundant:
                    isoform.isoformID = isoformID
                    isoformName = '%s_%d_%d' % (chrom,
                                                geneID,
                                                isoform.isoformID)
                    DNARecord = SeqRecord(isoform.dnaSeq,
                                            id=isoformName,
                                            description='''
                                            DNA sequence from predicted
                                            gene models of chickens line
                                            6 & 7
                                            ''')
                    isoformDNASeqs.append(DNARecord)

                    if isoform.frame:
                        proteinRecord = SeqRecord(isoform.proteinSeq,
                                                    id=isoformName,
                                                    description='''
                                                    protein sequence
                                                    from predicted
                                                    gene models of
                                                    chickens line 6 & 7
                                                    ''')
                        RNARecord = SeqRecord(isoform.mrnaSeq,
                                                id=isoformName,
                                                description='''
                                                mRNA sequence from predicted
                                                gene models of chickens line
                                                6 & 7
                                                ''')
                        isoformProteinSeqs.append(proteinRecord)
                        isoformRNASeqs.append(RNARecord)
                        isoformID += 1

                if n > 0 and n % 1000 == 0:
                    print >> sys.stderr, '...', n, 'transcripts done.'

    print >> sys.stderr, 'total genes = ', len(allGenes)
    print >> sys.stderr, 'Writing gene models to file...'
    writeBEDFile(allGenes, options.basename)
    print >> sys.stderr, 'Writing DNA sequences to file...'
    SeqIO.write(isoformDNASeqs, options.basename + '.dnas.fa', 'fasta')
    print >> sys.stderr, 'Writing RNA sequences to file...'
    SeqIO.write(isoformRNASeqs, options.basename + '.mrnas.fa', 'fasta')
    print >> sys.stderr, 'Writing protein sequences to file...'
    SeqIO.write(isoformProteinSeqs, options.basename + '.proteins.fa', 'fasta')

if __name__ == '__main__':
    (options, args) = parser.parse_args()
    main(options, args)
