import csv
import unittest
import bed2gff

TEST_FILE = 'COPS7A.models.bed'

class TestBed2Gff(unittest.TestCase):
    def setUp(self):
        self.allTranscripts = []
        for row in csv.reader(open(TEST_FILE), dialect='excel-tab'):
            transcript = bed2gff.parseBed(row)
            self.allTranscripts.append(transcript)

        '''Transcript Features
        
        1.Start Codon
        2.Stop Codon
        3.Strand
        4.Transcript name
        5. Chromosome
        6. Gene name
        7. Chromosome start
        8. Chrosome end
        9. Exons Starts
        10. Exon Sizes

        '''

    def testTranscriptFeatures(self):
        transcript = self.allTranscripts[0]
        self.assertEqual(transcript.chromosome, 'chr1')
        self.assertEqual(transcript.chromStart, 80303035)
        self.assertEqual(transcript.chromEnd, 80306166)
        self.assertEqual(transcript.strand, '+')
        self.assertEqual(transcript.transcriptName, '1')
        self.assertEqual(transcript.geneName, '8')
        self.assertEqual(transcript.exonStarts, [0, 352, 879, 1043, 1705, 2090, 2521])
        self.assertEqual(transcript.exonSizes, [260, 76, 89, 203, 106, 152, 610])
        #self.assertEqual(transcript.startCodon, (80303046, 80303048))
        #self.assertEqual(transcript.stopCodon, (80305594, 80305596))
        self.assertEqual((80303035, 80303045, 'UTR'), transcript.exons[0])
        self.assertEqual((80303046, 80303048, 'startCodon'), transcript.exons[1])
        #self.assertEqual((???,???, 'CDS'), transcript.exons[0])
        self.assertEqual(9, len(transcript.exons))

        transcript = self.allTranscripts[1]
        self.assertEqual(transcript.chromosome, 'chr1')
        self.assertEqual(transcript.chromStart, 80303035)
        self.assertEqual(transcript.chromEnd, 80306166)
        self.assertEqual(transcript.strand, '-')
        self.assertEqual(transcript.transcriptName, '1')
        self.assertEqual(transcript.geneName, '8')
        self.assertEqual(transcript.exonStarts, [0, 352, 879, 1043, 1705, 2090, 2521])
        self.assertEqual(transcript.exonSizes, [260, 76, 89, 203, 106, 152, 610])
        #self.assertEqual(transcript.startCodon, (80305594, 80305596))
        #self.assertEqual(transcript.stopCodon, (80303046, 80303048))
        self.assertEqual(9, len(transcript.exons))

        transcript = self.allTranscripts[2]
        self.assertEqual(transcript.chromosome, 'chr1')
        self.assertEqual(transcript.chromStart, 80303035)
        self.assertEqual(transcript.chromEnd, 80306166)
        self.assertEqual(transcript.strand, '.')
        self.assertEqual(transcript.transcriptName, '1')
        self.assertEqual(transcript.geneName, '8')
        self.assertEqual(transcript.exonStarts, [0, 352, 879, 1043, 1705, 2090, 2521])
        self.assertEqual(transcript.exonSizes, [260, 76, 89, 203, 106, 152, 610])

        self.assertEqual(7, len(transcript.exons))
        for start, end, seqType in transcript.exons:
            self.assertEqual('exon', seqType)
