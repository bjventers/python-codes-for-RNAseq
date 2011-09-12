import sys
import unittest

import bed2gff

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


class TestParseBed(unittest.TestCase):
    def setUp(self):
        self.row = ['chr1',
                    80303035,
                    80306166,
                    'chr1:8.1',
                    1000,
                    '+',
                    80303046,
                    80305596,
                    '0,0,0',
                    '7',
                    '260,76,89,203,106,152,610',
                    '0,352,879,1043,1705,2090,2521',
                    ]

    def testTranscriptFeaturesNormalCompletePositiveStrand(self):
        transcript = bed2gff.parseBed(self.row)

        self.assertEqual(transcript.chromosome, 'chr1')
        self.assertEqual(transcript.chromStart, 80303035)
        self.assertEqual(transcript.chromEnd, 80306166)
        self.assertEqual(transcript.strand, '+')
        self.assertEqual(transcript.transcriptName, '1')
        self.assertEqual(transcript.geneName, '8')
        self.assertEqual(transcript.exonStarts, [0, 352, 879, 1043, 1705, 2090, 2521])
        self.assertEqual(transcript.exonSizes, [260, 76, 89, 203, 106, 152, 610])


class TestConstructExons(unittest.TestCase):
    def setUp(self):
        self.row = ['chr1',
                    80303035,
                    80306166,
                    'chr1:8.1',
                    1000,
                    '+',
                    80303046,
                    80305596,
                    '0,0,0',
                    '7',
                    '260,76,89,203,106,152,610',
                    '0,352,879,1043,1705,2090,2521',
                    ]

    def testTranscriptFeaturesNormalCompletePositiveStrand(self):
        transcript = bed2gff.parseBed(self.row)

        #self.assertEqual(11, len(transcript.exons))

        self.assertEqual((80303035, 80303045, 'UTR'), transcript.exons[0])
        self.assertEqual((80303046, 80303048, 'startCodon'), transcript.exons[1])
        self.assertEqual((80305594, 80305596, 'stopCodon'), transcript.exons[-2])

        for start, end, seqType in transcript.exons[2:-2]:
            self.assertEqual(seqType, 'CDS')

    def testTranscriptFeaturesNormalCompleteNegativeStrand(self):
        self.row[5] = '-'  # change strand
        transcript = bed2gff.parseBed(self.row)

        #self.assertEqual(11, len(transcript.exons))

        self.assertEqual((80303035, 80303045, 'UTR'), transcript.exons[0])
        self.assertEqual((80303046, 80303048, 'stopCodon'), transcript.exons[1])
        self.assertEqual((80305594, 80305596, 'startCodon'), transcript.exons[-2])
        self.assertEqual((80305597, 80306166,'UTR'), transcript.exons[-1])

        for start, end, seqType in transcript.exons[2:-2]:
            self.assertEqual(seqType, 'CDS')

    def testTranscriptFeaturesNormalCompleteNoStrand(self):
        self.row[5] = '.'  # change strand
        transcript = bed2gff.parseBed(self.row)

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

    def testTranscriptFeaturesSplitStartCodonPositiveStrand(self):
        self.row[6] = 80303294  # change thickStart
        transcript = bed2gff.parseBed(self.row)

        self.assertEqual(11, len(transcript.exons))

        self.assertEqual((80303035, 80303293, 'UTR'), transcript.exons[0])
        self.assertEqual((80303294, 80303295, 'startCodon'), transcript.exons[1])
        self.assertEqual((80303387, 80303387, 'startCodon'), transcript.exons[2])
        self.assertEqual((80303388, 80303463, 'CDS'), transcript.exons[3])


        for start, end, seqType in transcript.exons[3:-2]:
            self.assertEqual(seqType, 'CDS')

        self.assertEqual('stopCodon', transcript.exons[-2][-1])
        self.assertEqual('UTR', transcript.exons[-1][-1])

    def testTranscriptFeaturesSplitStopCodonNegativeStrand(self):
        self.row[6] = 80303294  # modifythickStart
        self.row[5] = '-'  # modify strand
        transcript = bed2gff.parseBed(self.row)

        self.assertEqual(11, len(transcript.exons))

        self.assertEqual((80303035, 80303293, 'UTR'), transcript.exons[0])
        self.assertEqual((80303294, 80303295, 'stopCodon'), transcript.exons[1])
        self.assertEqual((80303387, 80303387, 'stopCodon'), transcript.exons[2])
        self.assertEqual((80303388, 80303463, 'CDS'), transcript.exons[3])


        for start, end, seqType in transcript.exons[3:-2]:
            self.assertEqual(seqType, 'CDS')

        self.assertEqual('startCodon', transcript.exons[-2][-1])
        self.assertEqual('UTR', transcript.exons[-1][-1])

    def testTranscriptFeaturesSplitStopCodonPositiveStrand(self):
        self.row[7] = 80305556
        transcript = bed2gff.parseBed(self.row)

        #self.assertEqual(11, len(transcript.exons))

        self.assertEqual((80303035, 80303045, 'UTR'), transcript.exons[0])
        self.assertEqual((80303046, 80303048, 'startCodon'), transcript.exons[1])
        self.assertEqual((80305276, 80305277, 'stopCodon'), transcript.exons[-3])
        self.assertEqual((80305556, 80305556, 'stopCodon'), transcript.exons[-2])
        self.assertEqual((80305557, 80306166, 'UTR'), transcript.exons[-1])

        for start, end, seqType in transcript.exons[2:-3]:
            self.assertEqual(seqType, 'CDS')

    def testTranscriptFeaturesSplitStartCodonNegativeStrand(self):
        self.row[7] = 80305556
        self.row[5] = '-'
        transcript = bed2gff.parseBed(self.row)

        #self.assertEqual(11, len(transcript.exons))

        self.assertEqual((80303035, 80303045, 'UTR'), transcript.exons[0])
        self.assertEqual((80303046, 80303048, 'stopCodon'), transcript.exons[1])
        self.assertEqual((80305276, 80305277, 'startCodon'), transcript.exons[-3])
        self.assertEqual((80305556, 80305556, 'startCodon'), transcript.exons[-2])
        self.assertEqual((80305557, 80306166, 'UTR'), transcript.exons[-1])

        for start, end, seqType in transcript.exons[2:-3]:
            self.assertEqual(seqType, 'CDS')

    def testTranscriptFeaturesSplitStartStopCodonPositiveStrand(self):
        self.row[6] = 80303294  # change thickStart
        self.row[7] = 80305556
        transcript = bed2gff.parseBed(self.row)

        self.assertEqual(11, len(transcript.exons))

        self.assertEqual((80303035, 80303293, 'UTR'), transcript.exons[0])
        self.assertEqual((80303294, 80303295, 'startCodon'), transcript.exons[1])
        self.assertEqual((80303387, 80303387, 'startCodon'), transcript.exons[2])
        self.assertEqual((80303388, 80303463, 'CDS'), transcript.exons[3])

        self.assertEqual((80305276, 80305277, 'stopCodon'), transcript.exons[-3])
        self.assertEqual((80305556, 80305556, 'stopCodon'), transcript.exons[-2])
        self.assertEqual((80305557, 80306166, 'UTR'), transcript.exons[-1])

        for start, end, seqType in transcript.exons[3:-3]:
            self.assertEqual(seqType, 'CDS')

    def testTranscriptFeaturesNormalPositiveStrandNo_5_UTR(self):
        self.row[6] = 80303035  # change thickStart
        transcript = bed2gff.parseBed(self.row)

        self.assertEqual(10, len(transcript.exons))

        self.assertEqual((80303035, 80303037, 'startCodon'), transcript.exons[0])

        for start, end, seqType in transcript.exons[1:-2]:
            self.assertEqual(seqType, 'CDS')

    def testTranscriptFeaturesNormalNegativeStrandNo_3_UTR(self):
        self.row[5] = '-'  # change strand
        self.row[6] = 80303035
        transcript = bed2gff.parseBed(self.row)

        #self.assertEqual(11, len(transcript.exons))

        self.assertEqual((80303035, 80303037, 'stopCodon'), transcript.exons[0])
        self.assertEqual((80305594, 80305596, 'startCodon'), transcript.exons[-2])
        self.assertEqual((80305597, 80306166,'UTR'), transcript.exons[-1])

        for start, end, seqType in transcript.exons[1:-2]:
            self.assertEqual(seqType, 'CDS')


class TestNewThickEnd(unittest.TestCase):
    def setUp(self):
        self.row = ['chr1',
                    80303035,
                    80306166,
                    'chr1:8.1',
                    1000,
                    '+',
                    80303046,
                    80305596,
                    '0,0,0',
                    '7',
                    '260,76,89,203,106,152,610',
                    '0,352,879,1043,1705,2090,2521',
                    ]

    def testNormalThickEnd(self):
        transcript = bed2gff.parseBed(self.row)
        self.assertEqual(80305594, transcript._newThickEnd())

    def testSplitThickEnd1(self):
        self.row[7] = 80305556
        transcript = bed2gff.parseBed(self.row)
        self.assertEqual(80305276, transcript._newThickEnd())

    def testSplitThickEnd2(self):
        self.row[7] = 80305557
        transcript = bed2gff.parseBed(self.row)
        self.assertEqual(80305277, transcript._newThickEnd())

    def testSplitThickEnd3(self):
        self.row[7] = 80305277
        transcript = bed2gff.parseBed(self.row)
        self.assertEqual(80305275, transcript._newThickEnd())

    def testSplitThickEnd4(self):
        self.row[7] = 80305558
        transcript = bed2gff.parseBed(self.row)
        self.assertEqual(80305556, transcript._newThickEnd())
