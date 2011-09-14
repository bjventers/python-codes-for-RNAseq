'''
    This script converts BED format to GFF format.
    Note: the first base in BED format is numbered 0.
    Author: Likit Preeyanon
'''

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
        self.thickStart = int(row[6]) + 1
        self.thickEnd = int(row[7])
        self.exons = self._constructExons()
        self.frames = self._getFrame()

    def _newThickEnd(self):
        '''Return thickEnd - 2 as a new thickEnd.
        
        This function shift thickEnd to prior exon if needed.

        '''

        for i in range(len(self.exonStarts)):
            start = self.exonStarts[i] + self.chromStart + 1
            end = self.exonSizes[i] + start - 1
            if (start <= self.thickEnd) and (end >= self.thickEnd):
                if self.thickEnd - 2 < start:
                    theRest = 3 - ((self.thickEnd - start) + 1)
                    try:
                        prevStart = self.exonStarts[i - 1] + self.chromStart + 1
                        prevEnd = self.exonSizes[i - 1] + prevStart - 1
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
        start = self.exonStarts[i] + self.chromStart + 1
        end = self.exonSizes[i] + start - 1
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
                            UTRType = '5UTR' if self.strand == '+' else '3UTR'
                            exons.append((start, newEnd, UTRType))

                    if (end - thickStart) + 1 >= remainSSCodonBase:
                        # number of bases in a start/stop codon in this exon.
                        SSCodonBase = remainSSCodonBase
                    else:
                        SSCodonBase = remainSSCodonBase - (end - thickStart)

                    remainSSCodonBase = remainSSCodonBase - SSCodonBase

                    if self.strand == '+':
                        exons.append((thickStart,
                                    thickStart + SSCodonBase - 1,
                                    'start_codon'))
                    else:
                        exons.append((thickStart,
                                    thickStart + SSCodonBase - 1,
                                    'stop_codon'))

                    if self.strand == '+':
                        remainBase = (end + 1) - (thickStart)
                    else:
                        remainBase = (end + 1) - (thickStart + SSCodonBase)

                    if remainBase > 0:
                        '''
                            If there are some bases after a start/stop codon,
                            assign them to CDS.
                        '''
                        if self.strand == '+':
                            exons.append((thickStart, end, 'CDS'))
                        else:
                            exons.append((thickStart + SSCodonBase, end, 'CDS'))
                else:
                    exons.append((start, end, 'CDS'))

            start = self.exonStarts[i] + self.chromStart + 1
            end = self.exonSizes[i] + start - 1

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

                        if self.strand == '+':
                            newEnd = thickEnd - 1
                        else:
                            newEnd = thickEnd + 2

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
                                    'stop_codon'))
                    else:
                        if (start, thickEnd + SSCodonBase - 1, 'CDS') not \
                                in exons:
                            exons.append((start,
                                            thickEnd + SSCodonBase - 1,
                                            'CDS'))
                        exons.append((thickEnd,
                                        thickEnd + SSCodonBase - 1,
                                        'start_codon'))

                    remainBase = (end + 1) - (thickEnd + SSCodonBase)

                    if remainBase > 0:
                        '''
                            If there are some bases after a start/stop codon,
                            assign them to UTR.
                        '''
                        UTRType = '5UTR' if self.strand == '-' else '3UTR'
                        exons.append((thickEnd + SSCodonBase, end, UTRType))
                else:
                    UTRType = '5UTR' if self.strand == '-' else '3UTR'
                    exons.append((start, end, UTRType))

            try:
                start = self.exonStarts[i] + self.chromStart + 1
                end = self.exonSizes[i] + start - 1
            except IndexError:
                break

        return exons

    def _getFrame(self):

        frames = []
        frame = 0
        prevExonSize = 0

        for i in range(len(self.exons)):
            if self.exons[i][2] == '5UTR' or self.exons[i][2] == '3UTR':
                frames.append(None)
            else:
                frame = (3 - ((prevExonSize - frame) % 3)) % 3
                frames.append(frame)
                prevExonSize = self.exons[i][1] - self.exons[i][0] + 1

        return frames

def parseBed(row):
    assert len(row) == 12, "Number of row must be 12"
    return Transcript(row)

def printGTF(transcript):
    for i in range(len(transcript.exons)):
        start, end, seqType = transcript.exons[i]

        if transcript.frames[i] != None:
            frame = str(transcript.frames[i])
        else:
            frame = '.'

        print '%s\tRNA-Seq\t%s\t%d\t%d\t.\t%s\t%s\t \
                gene_name "%s"; transcript_name "%s";exon_number "%d"' \
                                            % (transcript.chromosome,
                                                seqType,
                                                start,
                                                end,
                                                transcript.strand,
                                                frame,
                                                transcript.geneName,
                                                transcript.transcriptName,
                                                i + 1)


if __name__=='__main__':

    fp = open(sys.argv[1])

    for line in csv.reader(fp, dialect='excel-tab'):
        transcript = parseBed(line)
        printGTF(transcript)
