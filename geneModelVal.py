import sys
import csv
from pygr import seqdb
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC

class Isoform(object):
    """docstring for Isoform"""
    def __init__(self, genome, *args):
        self.chrom = args[0]
        self.chromStart = int(args[1])
        self.chromEnd = int(args[2])
        self.name = args[3]
        self.score = args[4]
        self.strand = args[5]
        self.startCodon = int(args[6])
        self.stopCodon = int(args[7])
        self.itemRgb = args[8]
        self.blockCount = args[9]
        self.blockSizes = [int(s) for s in args[10].split(',')]
        self.blockStarts = [int(s) for s in args[11].split(',')]
        self.exons = self.getExons()
        self.dnaSeq = self.getDnaSequence()
        self.proteinSeq = self.getProteinSequence()

    def getDnaSequence(self):
        seq = ''
        for start, stop in self.exons:
            seq += str(genome[self.chrom][start:stop])
        bioSeq = Seq(seq.upper(), IUPAC.unambiguous_dna)
        if self.strand == '+':
            return bioSeq 
        else:
            return bioSeq.reverse_complement()

    def getProteinSequence(self):
        orfStart = self.startCodon-self.chromStart
        orfEnd = self.stopCodon-self.chromStart
        self.orf = self.dnaSeq[orfStart:orfEnd]
        print 'ORF start, end = ', orfStart, orfEnd
        return self.orf.translate()

    def getExons(self):
        exons = []
        for i in range(len(self.blockStarts)):
            start = self.blockStarts[i]
            stop = self.blockSizes[i] + self.blockStarts[i]
            exons.append((start,stop))
        return exons
        
allGenes = {}

genome = seqdb.SequenceFileDB(sys.argv[2])
modelReader = csv.reader(open(sys.argv[1], 'rb'), delimiter='\t')
for row in modelReader:
    name = row[3]
    chrom = row[0]
    chrom, geneID, isoID, _, frame = name.split('_')
    isoform = Isoform(genome, *row)
    if chrom not in allGenes:
        allGenes[chrom] = {}
        allGenes[chrom][geneID] = [isoform]
    else:
        try:
            allGenes[chrom][geneID].append(isoform)
        except KeyError:
            allGenes[chrom][geneID] = [isoform]

print 'total gene:', len(allGenes)
print isoform.name, isoform.chromStart, isoform.chromEnd, isoform.chrom, isoform.startCodon, isoform.stopCodon
print isoform.exons, isoform.strand
print isoform.orf
print isoform.proteinSeq
