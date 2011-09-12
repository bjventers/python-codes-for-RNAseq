import sys
import csv

class Transcript(object):
    def __init__(self, row):
        self.chromStart = int(row[1])
        self.chromEnd = int(row[2])
        self.chromosome = row[3].split(':')[0]  
        self.geneName, self.transcriptName = row[3].split(':')[1].split('.')
        self.strand = row[5]
        self.exonStarts = [int(i) for i in row[-1].split(',')]
        self.exonSizes = [int(i) for i in row[-2].split(',')]
        self.thickStart = int(row[6])
        self.thickEnd = int(row[7])
        self.exons = self._constructExons()

    def _newThickEnd(self):

        for i in range(len(self.exonStarts)):
            start = self.exonStarts[i] + self.chromStart
            end = self.exonSizes[i] + start
            if (start <= self.thickEnd) and (end >= self.thickEnd):
                if self.thickEnd - 2 < start:
                    theRest = 3 - ((self.thickEnd - start) + 1)
                    try:
                        prevStart = self.exonStarts[i - 1] + self.chromStart
                        prevEnd = self.exonSizes[i - 1] + prevStart
                    except IndexError:
                        return None
                    else:
                        return prevEnd - theRest + 1
                else:
                    return self.thickEnd - 2


    def _constructExons(self):
        thickStart = self.thickStart
        thickEnd = self._newThickEnd()

        exons = []

        '''Number of base needed to complete a start/stop codon.'''
        remainSSCodonBase = 3

        i = 0
        start = self.exonStarts[i] + self.chromStart
        end = self.exonSizes[i] + start
        while (end < thickEnd):
            i += 1
            if remainSSCodonBase > 0 and remainSSCodonBase < 3:
                '''
                    If a start/stop codon is found but not complete,
                    assign a number of bases in the next exon to
                    complete a start/stop codon.
                '''
                thickStart = start

            if self.strand == '.':
                exons.append((start, end, 'exon'))
            else:
                if (start <= thickStart) and (end >= thickStart):
                    '''
                        If the exon contains start or stop codon,
                        all bases before a start/stop codon belong
                        to UTR and all bases after a start/stop codon
                        belong to CDS.
                    '''


                    if remainSSCodonBase == 3:
                        if start <= thickStart - 1:
                            newEnd = thickStart - 1
                            exons.append((start, newEnd, 'UTR'))

                    if (end - thickStart) + 1 >= remainSSCodonBase:
                        # number of bases in a start/stop codon in this exon.
                        SSCodonBase = remainSSCodonBase
                    else:
                        SSCodonBase = remainSSCodonBase - (end - thickStart)

                    remainSSCodonBase = remainSSCodonBase - SSCodonBase

                    if self.strand == '+':
                        exons.append((thickStart,
                                    thickStart + SSCodonBase - 1,
                                    'startCodon'))
                    else:
                        exons.append((thickStart,
                                    thickStart + SSCodonBase - 1,
                                    'stopCodon'))

                    remainBase = (end + 1) - (thickStart + SSCodonBase)

                    if remainBase > 0:
                        '''
                            If there are some bases after a start/stop codon,
                            assign them to CDS.
                        '''
                        exons.append((thickStart + SSCodonBase, end, 'CDS'))
                else:
                    exons.append((start, end, 'CDS'))

            start = self.exonStarts[i] + self.chromStart
            end = self.exonSizes[i] + start

        '''Number of base needed to complete a start/stop codon.'''
        remainSSCodonBase = 3

        while i < len(self.exonStarts):
            i += 1
            if remainSSCodonBase > 0 and remainSSCodonBase < 3:
                '''
                    If a start/stop codon is found but not complete,
                    assign a number of bases in the next exon to
                    complete a start/stop codon.
                '''
                thickEnd = start

            if self.strand == '.':
                exons.append((start, end, 'exon'))
            else:
                if (start <= thickEnd) and (end >= thickEnd):
                    '''
                        If the exon contains start or stop codon,
                        all bases before a start/stop codon belong
                        to UTR and all bases after a start/stop codon
                        belong to CDS.
                    '''
                    if remainSSCodonBase == 3:

                        newEnd = thickEnd - 1
                        exons.append((start, newEnd, 'CDS'))

                    if (end - thickEnd) + 1 >= remainSSCodonBase:
                        # number of bases in a start/stop codon in this exon.
                        SSCodonBase = remainSSCodonBase
                    else:
                        SSCodonBase = remainSSCodonBase - (end - thickEnd)

                    remainSSCodonBase = remainSSCodonBase - SSCodonBase

                    if self.strand == '+':
                        exons.append((thickEnd,
                                    thickEnd + SSCodonBase - 1,
                                    'stopCodon'))
                    else:
                        exons.append((thickEnd,
                                    thickEnd + SSCodonBase - 1,
                                    'startCodon'))

                    remainBase = (end + 1) - (thickEnd + SSCodonBase)

                    if remainBase > 0:
                        '''
                            If there are some bases after a start/stop codon,
                            assign them to CDS.
                        '''
                        exons.append((thickEnd + SSCodonBase, end, 'UTR'))
                else:
                    exons.append((start, end, 'UTR'))

            try:
                start = self.exonStarts[i] + self.chromStart
                end = self.exonSizes[i] + start
            except IndexError:
                break

        return exons


def parseBed(row):
    assert len(row) == 12, "Number of row must be 12"
    transcript = Transcript(row)
    return transcript
