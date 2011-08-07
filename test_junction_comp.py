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

    def test_find_match_matched(self):
        key = 'chr1:154000-230000'
        mocker = Mocker()

        junction = mocker.mock()
        junction.coverage
        mocker.result(40)
        mocker.count(1, None)

        container = mocker.mock()
        container.keys()
        mocker.result([key])

        container[key]
        mocker.result(junction)
        mocker.count(1, None)

        mocker.replay()

        self.common, self.diff = jc.findMatch(container, container)
        self.assertEqual(self.common.keys(), ['chr1:154000-230000'])
        self.assertEqual(self.diff.keys(), [])

        mocker.restore()
        mocker.verify()

    def test_find_match_not_matched(self):
        key1 = 'chr1:154000-230000'
        mocker = Mocker()

        junction = mocker.mock()
        junction.coverage
        mocker.result(40)

        container1 = mocker.mock()
        container2 = mocker.mock()

        container1.keys()
        mocker.result([key1])

        container1[key1]
        mocker.result(junction)

        container2[key1]
        mocker.throw(KeyError)
        mocker.count(1)

        mocker.replay()

        self.common, self.diff = jc.findMatch(container1, container2)
        self.assertEqual(self.common.keys(), [])
        self.assertEqual(self.diff.keys(), [key1])

        mocker.restore()
        mocker.verify()

    def test_find_match_match_not_match(self):
        key1 = ['chr1:154000-230000', 'chr1:155000-230000']
        key2 = ['chr1:154000-230000']
        mocker = Mocker()
        junction = mocker.mock()
        junction.coverage
        mocker.result(40)
        mocker.count(1, None)

        container1 = mocker.mock()
        container2 = mocker.mock()

        container1.keys()
        mocker.result(key1)

        container1[key1[0]]
        mocker.result(junction)
        mocker.count(1, None)

        container1[key1[1]]
        mocker.result(junction)
        mocker.count(1, None)

        container2[key1[0]]
        mocker.result(junction)
        mocker.count(1, None)

        container2[key1[1]]
        mocker.throw(KeyError)
        mocker.count(1)

        mocker.replay()
        self.common, self.diff = jc.findMatch(container1, container2)
        self.assertEqual(self.common.keys(), [key1[0]])
        self.assertEqual(self.diff.keys(), [key1[1]])

        mocker.restore()
        mocker.verify()


class TestJunctionScan(TestCase):

    def test_shared_exon(self):
        mocker = Mocker()
        self.model = mocker.mock()
        self.model.blockStarts
        mocker.result([100, 400, 700])
        self.model.blockEnds
        mocker.result([200, 500, 800])
        self.model.chromName
        mocker.result('chr1')
        mocker.replay()
        self.junctions1 = {200: 400, 500: 700}
        self.junctions2 = {200: 700}
        self.sharedExons = jc.scanJunctions(self.model,
                                self.junctions1, self.junctions2)

        for j in self.sharedExons:
            self.assertEqual(j, [400, 720])


if __name__ == '__main__':
    main()
