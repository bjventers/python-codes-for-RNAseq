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

        self.assertEqual(11, len(transcript.exons))

        print >> sys.stderr, '\n'
        for i in range(len(transcript.exons)):
            print >> sys.stderr, i + 1, transcript.exons[i]

        self.assertEqual((80303036, 80303046, '5UTR'), transcript.exons[0])
        self.assertEqual((80303047, 80303049, 'start_codon'), transcript.exons[1])
        self.assertEqual((80303047, 80303295, 'CDS'), transcript.exons[2])
        self.assertEqual((80305594, 80305596, 'stop_codon'), transcript.exons[-2])

        for start, end, seqType in transcript.exons[2:-2]:
            self.assertEqual(seqType, 'CDS')

    def testTranscriptFeaturesNormalCompleteNegativeStrand(self):
        self.row[5] = '-'  # change strand
        transcript = bed2gff.parseBed(self.row)

        self.assertEqual(11, len(transcript.exons))

        self.assertEqual((80303036, 80303046, '3UTR'), transcript.exons[0])
        self.assertEqual((80303047, 80303049, 'stop_codon'), transcript.exons[1])
        self.assertEqual((80303050, 80303295, 'CDS'), transcript.exons[2])
        self.assertEqual((80305557, 80305596, 'CDS'), transcript.exons[-3])
        self.assertEqual((80305594, 80305596, 'start_codon'), transcript.exons[-2])
        self.assertEqual((80305597, 80306166,'5UTR'), transcript.exons[-1])

        for start, end, seqType in transcript.exons[2:-2]:
            self.assertEqual(seqType, 'CDS')

    def testTranscriptFeaturesNormalCompleteNoStrand(self):
        self.row[5] = '.'  # change strand
        transcript = bed2gff.parseBed(self.row)

        self.assertEqual(7, len(transcript.exons))

        for start, end, seqType in transcript.exons:
            self.assertEqual('exon', seqType)

    def testTranscriptFeaturesSplitStartCodonPositiveStrand(self):
        self.row[6] = 80303293  # change thickStart
        transcript = bed2gff.parseBed(self.row)

        self.assertEqual(12, len(transcript.exons))

        self.assertEqual((80303036, 80303293, '5UTR'), transcript.exons[0])
        self.assertEqual((80303294, 80303295, 'start_codon'), transcript.exons[1])
        self.assertEqual((80303294, 80303295, 'CDS'), transcript.exons[2])
        self.assertEqual((80303388, 80303388, 'start_codon'), transcript.exons[3])
        self.assertEqual((80303388, 80303463, 'CDS'), transcript.exons[4])


        for start, end, seqType in transcript.exons[4:-2]:
            self.assertEqual(seqType, 'CDS')

        self.assertEqual('stop_codon', transcript.exons[-2][-1])
        self.assertEqual('3UTR', transcript.exons[-1][-1])

    def testTranscriptFeaturesSplitStopCodonNegativeStrand(self):
        self.row[6] = 80303293  # modifythickStart
        self.row[5] = '-'  # modify strand
        transcript = bed2gff.parseBed(self.row)

        self.assertEqual(11, len(transcript.exons))

        self.assertEqual((80303036, 80303293, '3UTR'), transcript.exons[0])
        self.assertEqual((80303294, 80303295, 'stop_codon'), transcript.exons[1])
        self.assertEqual((80303388, 80303388, 'stop_codon'), transcript.exons[2])
        self.assertEqual((80303389, 80303463, 'CDS'), transcript.exons[3])


        for start, end, seqType in transcript.exons[3:-2]:
            self.assertEqual(seqType, 'CDS')

        self.assertEqual('start_codon', transcript.exons[-2][-1])
        self.assertEqual('5UTR', transcript.exons[-1][-1])

    def testTranscriptFeaturesSplitStopCodonPositiveStrand(self):
        self.row[7] = 80305557
        transcript = bed2gff.parseBed(self.row)

        self.assertEqual(11, len(transcript.exons))

        self.assertEqual((80303036, 80303046, '5UTR'), transcript.exons[0])
        self.assertEqual((80303047, 80303049, 'start_codon'), transcript.exons[1])
        self.assertEqual((80305276, 80305277, 'stop_codon'), transcript.exons[-3])
        self.assertEqual((80305557, 80305557, 'stop_codon'), transcript.exons[-2])
        self.assertEqual((80305558, 80306166, '3UTR'), transcript.exons[-1])

        for start, end, seqType in transcript.exons[2:-3]:
            self.assertEqual(seqType, 'CDS')

    def testTranscriptFeaturesSplitStartCodonNegativeStrand(self):
        self.row[7] = 80305557
        self.row[5] = '-'
        transcript = bed2gff.parseBed(self.row)

        self.assertEqual(13, len(transcript.exons))

        self.assertEqual((80303036, 80303046, '3UTR'), transcript.exons[0])
        self.assertEqual((80303047, 80303049, 'stop_codon'), transcript.exons[1])
        self.assertEqual((80305558, 80306166, '5UTR'), transcript.exons[-1])
        self.assertEqual((80305126, 80305277, 'CDS'), transcript.exons[-5])
        self.assertEqual((80305276, 80305277, 'start_codon'), transcript.exons[-4])
        self.assertEqual((80305557, 80305557, 'CDS'), transcript.exons[-3])
        self.assertEqual((80305557, 80305557, 'start_codon'), transcript.exons[-2])

        for start, end, seqType in transcript.exons[2:-5]:
            self.assertEqual(seqType, 'CDS')

    def testTranscriptFeaturesSplitStartStopCodonPositiveStrand(self):
        self.row[6] = 80303293  # change thickStart
        self.row[7] = 80305557
        transcript = bed2gff.parseBed(self.row)

        self.assertEqual(12, len(transcript.exons))

        self.assertEqual((80303036, 80303293, '5UTR'), transcript.exons[0])
        self.assertEqual((80303294, 80303295, 'start_codon'), transcript.exons[1])
        self.assertEqual((80303294, 80303295, 'CDS'), transcript.exons[2])
        self.assertEqual((80303388, 80303388, 'start_codon'), transcript.exons[3])
        self.assertEqual((80303388, 80303463, 'CDS'), transcript.exons[4])

        self.assertEqual((80305276, 80305277, 'stop_codon'), transcript.exons[-3])
        self.assertEqual((80305557, 80305557, 'stop_codon'), transcript.exons[-2])
        self.assertEqual((80305558, 80306166, '3UTR'), transcript.exons[-1])

        for start, end, seqType in transcript.exons[4:-3]:
            self.assertEqual(seqType, 'CDS')

    def testTranscriptFeaturesNormalPositiveStrandNo_5_UTR(self):
        self.row[6] = 80303035  # change thickStart
        transcript = bed2gff.parseBed(self.row)

        self.assertEqual(10, len(transcript.exons))

        self.assertEqual((80303036, 80303038, 'start_codon'), transcript.exons[0])

        for start, end, seqType in transcript.exons[1:-2]:
            self.assertEqual(seqType, 'CDS')

    def testTranscriptFeaturesNormalPositiveStrandNo_3_UTR(self):
        self.row[7] = 80306166  # change thickStart
        transcript = bed2gff.parseBed(self.row)

        self.assertEqual(10, len(transcript.exons))

        self.assertEqual((80306164, 80306166, 'stop_codon'), transcript.exons[-1])

        for start, end, seqType in transcript.exons[2:-1]:
            self.assertEqual(seqType, 'CDS')

    def testTranscriptFeaturesNormalNegativeStrandNo_3_UTR(self):
        self.row[5] = '-'  # change strand
        self.row[6] = 80303035
        transcript = bed2gff.parseBed(self.row)

        self.assertEqual(10, len(transcript.exons))

        self.assertEqual((80303036, 80303038, 'stop_codon'), transcript.exons[0])
        self.assertEqual((80305594, 80305596, 'start_codon'), transcript.exons[-2])
        self.assertEqual((80305597, 80306166,'5UTR'), transcript.exons[-1])

        for start, end, seqType in transcript.exons[1:-2]:
            self.assertEqual(seqType, 'CDS')

    def testTranscriptFeaturesNormalNegativeStrandNo_5_UTR(self):
        self.row[5] = '-'  # change strand
        self.row[7] = 80306166
        transcript = bed2gff.parseBed(self.row)

        self.assertEqual(10, len(transcript.exons))

        self.assertEqual((80306164, 80306166, 'start_codon'), transcript.exons[-1])

        for start, end, seqType in transcript.exons[2:-1]:
            self.assertEqual(seqType, 'CDS')

    def testTranscriptFeaturesNormalPositiveStrandNo_5_3_UTR(self):
        self.row[6] = 80303035  # change thickStart
        self.row[7] = 80306166  # change thickEnd
        transcript = bed2gff.parseBed(self.row)

        self.assertEqual(9, len(transcript.exons))

        self.assertEqual((80303036, 80303038, 'start_codon'), transcript.exons[0])
        self.assertEqual((80306164, 80306166, 'stop_codon'), transcript.exons[-1])

        for start, end, seqType in transcript.exons[1:-1]:
            self.assertEqual(seqType, 'CDS')

    def testTranscriptFeaturesNormalNegativeStrandNo_5_3_UTR(self):
        self.row[5] = '-'  # change strand
        self.row[6] = 80303035  # change thickStart
        self.row[7] = 80306166  # change thickEnd

        transcript = bed2gff.parseBed(self.row)

        self.assertEqual(9, len(transcript.exons))

        self.assertEqual((80303036, 80303038, 'stop_codon'), transcript.exons[0])
        self.assertEqual((80306164, 80306166, 'start_codon'), transcript.exons[-1])

        for start, end, seqType in transcript.exons[1:-1]:
            self.assertEqual(seqType, 'CDS')

    def testTranscriptFeaturesNormalPositiveStrandShort_5_3_UTR(self):
        self.row[6] = 80303036  # change thickStart
        self.row[7] = 80306165  # change thickEnd

        transcript = bed2gff.parseBed(self.row)

        '''
        print >> sys.stderr, '\n'
        for i in range(len(transcript.exons)):
            print >> sys.stderr, i + 1, transcript.exons[i]
        '''

        self.assertEqual(11, len(transcript.exons))

        self.assertEqual((80303036, 80303036,'5UTR'), transcript.exons[0])
        self.assertEqual((80306166, 80306166,'3UTR'), transcript.exons[-1])

    def testTranscriptFeaturesNormalNegativeStrandShort_5_3_UTR(self):
        self.row[5] = '-'
        self.row[6] = 80303036  # change thickStart
        self.row[7] = 80306165  # change thickEnd

        transcript = bed2gff.parseBed(self.row)

        self.assertEqual(11, len(transcript.exons))

        self.assertEqual((80303036, 80303036,'3UTR'), transcript.exons[0])
        self.assertEqual((80306166, 80306166,'5UTR'), transcript.exons[-1])


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
        self.row[7] = 80305557
        transcript = bed2gff.parseBed(self.row)
        self.assertEqual(80305276, transcript._newThickEnd())

    def testSplitThickEnd2(self):
        self.row[7] = 80305558
        transcript = bed2gff.parseBed(self.row)
        self.assertEqual(80305277, transcript._newThickEnd())

    def testSplitThickEnd3(self):
        self.row[7] = 80305277
        transcript = bed2gff.parseBed(self.row)
        self.assertEqual(80305275, transcript._newThickEnd())

    def testSplitThickEnd4(self):
        self.row[7] = 80305559
        transcript = bed2gff.parseBed(self.row)
        self.assertEqual(80305557, transcript._newThickEnd())


class TestGetFrame(unittest.TestCase):
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

    def testNormalCompletePositiveStrand(self):

        transcript = bed2gff.Transcript(self.row)

        frame = [None, 0, 0, 0, 2, 0, 1, 0, 1, 0, None]
        self.assertEqual(transcript.frames, frame)

    def testNormalCompleteNegativeStrand(self):
        self.row[5] = '-'
        transcript = bed2gff.Transcript(self.row)

        for i in range(len(transcript.frames)):
            print >> sys.stderr, transcript.exons[i], transcript.frames[i]

        frame = [None, 0, 0, 0, 2, 0, 1, 0, 1, 0, None]
        self.assertEqual(transcript.frames, frame)
