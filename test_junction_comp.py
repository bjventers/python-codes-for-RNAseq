from unittest import TestCase, main
import junction_comp as jc
from mocker import Mocker


class TestJunctionClass(TestCase):

    def test_construct(self):
        self.junc = jc.Junction('chr1', 154000, 230000,
                                230, '+', 'JUNC00001')
        self.assertEqual(self.junc.chrom, 'chr1')
        self.assertEqual(self.junc.start, 154000)
        self.assertEqual(self.junc.end, 230000)
        self.assertEqual(self.junc.coverage, 230)
        self.assertEqual(self.junc.name, 'JUNC00001')
        self.assertEqual(self.junc.strand, '+')

    def test_get_coord(self):
        mocker = Mocker()
        junc = mocker.mock()
        junc.chrom
        mocker.result('chr1')
        junc.start
        mocker.result(154000)
        junc.end
        mocker.result(230000)
        mocker.replay()
        self.assertEqual(jc.Junction.getCoord.im_func(junc),
                        'chr1:154000-230000')

    def test_str(self):
        mocker = Mocker()
        junc = mocker.mock()
        junc.getCoord()
        mocker.result('chr1:154000-230000')
        junc.name
        mocker.result('JUNC00001')
        mocker.replay()
        self.assertEqual(jc.Junction.__str__.im_func(junc),
                        'chr1:154000-230000, JUNC00001')


class TestGetJunction(TestCase):
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

    def test_get_junction(self):
        for junction in jc.getJunction(self.row):
            self.assertEqual(junction.start, 6477)
            self.assertEqual(junction.end, 9995)

        junctions = [junc for junc in jc.getJunction(self.row)]
        self.assertEqual(len(junctions), 1)

class TestFindMatch(TestCase):

    def test_find_match(self):
        key = 'chr1:154000-230000'
        mocker = Mocker()
        junction = mocker.mock()
        container = mocker.mock()
        container.keys()
        mocker.result([key])
        container[key]
        mocker.result(junction)
        mocker.replay()
        self.common = jc.findMatch(container, container)
        self.assertEqual(self.common.keys(), ['chr1:154000-230000'])


if __name__ == '__main__':
    main()
