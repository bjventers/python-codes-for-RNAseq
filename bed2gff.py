import sys
import csv

class Isoform(object):
    def __init__(self, *row):
        self.txStart = int(row[1])
        self.txEnd = int(row[2])
        self.chrom, self.geneID, self.isoformID = row[3].split('_')
        self.score = row[4]
        self.strand = row[5]

        if self.strand == '+':
            self.cdsStart = int(row[6])
            self.cdsEnd = int(row[7])
        elif self.strand == '-':
            self.cdsEnd = int(row[6])
            self.cdsStart = int(row[7])
        else:
            self.cdsStart = None
            self.cdsEnd = None

        self.itemRGB = row[8]
        self.exonCount = int(row[9])
        self.exonSizes = row[10].split(',')
        self.exonStarts = row[11].split(',')
        self.__getExonSets()

    def __getExonSets(self):
        self.exonSets = []
        for i in range(len(self.exonStarts)):
            start = int(self.exonStarts[i]) + self.txStart
            end = start + int(self.exonSizes[i])
            self.exonSets.append((start, end))
        
for row in csv.reader(open(sys.argv[1]), dialect='excel-tab'):
    isoform = Isoform(*row)
    break

print isoform.chrom, isoform.txStart, isoform.txEnd, isoform.strand
print isoform.exonSets
