import unittest
import sys

from pygr import seqdb
import genebuilder


class TestGetDnaSeq(unittest.TestCase):
    def setUp(self):

        genomeFile = '/Users/Likit/projects/mdv/data/chick.fa'
        self.genome = seqdb.SequenceFileDB(genomeFile, verbose=False)

        exons = [('chr1', 51035309, 51035430), ('chr1', 51062489, 51062516)]

        self.isoform = genebuilder.Isoform('chr1', '1', '0', exons, self.genome)

    def testGetDnaSeq(self):
        seq = '''GGTGACTCCACCACCTCCCTGGGCAGCCCATTCCAGTGCCTGACCATCCCTTTCAGAGAAACTTTTCCTAACACCCAACCTGAATCTTCCCTGCCGCGAAACTGAGGCCATCCCCTCTAGTCCTACTGCTAGTTATGTGGGAAAAGAG'''

        dnaSeq = self.isoform._getDnaSeq(self.genome)
        self.assertEqual(seq, str(dnaSeq.reverse_complement()))


    def testFindORF(self):

        genebuilder.findORF(self.isoform)

        mrnaSeq = 'GGUGACUCCACCACCUCCCUGGGCAGCCCAUUCCAGUGCCUGACCAUCCCUUUCAGAGAAACUUUUCCUAACACCCAACCUGAAUCUUCCCUGCCGCGAAACUGAGGCCAUCCCCUCUAGUCCUACUGCUAGUUAUGUGGGAAAAGAG'

        self.assertLess(self.isoform.frame, 0)
        self.assertEqual('MGSPFQCLTIPFRETFPNTQPESSLPRN', str(self.isoform.proteinSeq))
        self.assertEqual(mrnaSeq, str(self.isoform.mrnaSeq))
