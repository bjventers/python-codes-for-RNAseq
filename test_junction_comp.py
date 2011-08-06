from unittest import TestCase, main
import junction_comp as jc
from mocker import Mocker


class TestJunctionClass(TestCase):

    def setUp(self):
        self.junc = jc.Junction('chr1', 154000, 230000,
                                230, '+', 'JUNC00001')

    def test_construct(self):
        self.assertEqual(self.junc.chrom, 'chr1')
        self.assertEqual(self.junc.start, 154000)
        self.assertEqual(self.junc.end, 230000)
        self.assertEqual(self.junc.coverage, 230)
        self.assertEqual(self.junc.name, 'JUNC00001')
        self.assertEqual(self.junc.strand, '+')

    def test_get_coord(self):
        self.assertEqual(self.junc.getCoord(),
                        'chr1:154000-230000')

    def test_str(self):
        mocker = Mocker()
        junction = mocker.mock()
        junction.getCoord()
        mocker.result('chr1:154000-230000')
        mocker.replay()
        self.junc.getCoord = junction.getCoord
        self.assertEqual(str(self.junc),
                        'chr1:154000-230000, JUNC00001')


class TestGetJunctionStartEnd(TestCase):
    def setUp(self):
        self.row = ['chr1',     # chromosome name
                    '6419',     # chromosome start
                    '10054',    # chromosome end
                    'JUNC00080661',  # junction name
                    '5',        # coverage
                    '+',        # strand
                    '6419',     # thick start
                    '10054',    # thick end
                    '255,0,0',  # RGB
                    '2',        # block count
                    '58,59',    # block sizes
                    '0,3576',   # block starts
                    ]

    def test_getJunctionStartEnd(self):
        self.juncStart, self.juncEnd = jc.getJunctionStartEnd(self.row)
        self.assertEqual(self.juncStart, 6477)
        self.assertEqual(self.juncEnd, 9995)


if __name__ == '__main__':
    main()
