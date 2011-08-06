from unittest import TestCase
import junction_comp as jc

class testJunctionClass(TestCase):

    def test_construct(self):
        junc = jc.Junction('chr1', 154000, 230000, 230, 'JUNC00001')
        self.assertEqual(junc.chrom, 'chr1')
        self.assertEqual(junc.start, 154000)
        self.assertEqual(junc.end, 230000)
        self.assertEqual(junc.coverage, 230)
        self.assertEqual(junc.name, 'JUNC00001')

    def test_get_coord(self):
        self.assertEqual(junc.getCoord(), 'chr1:154000-230000')

    def test_str(self):
        self.assertEqual(str(junc), 'chr1:154000-230000, JUNC00001')
