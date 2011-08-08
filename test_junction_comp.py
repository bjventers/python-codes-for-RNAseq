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
        mocker.count(2, None)

        self.model.blockSizes
        mocker.result([100, 100, 100])
        mocker.count(1, None)

        self.model.chrom
        mocker.result('chr1')
        mocker.count(2, None)

        self.model.start
        mocker.result(0)
        mocker.count(1, None)

        self.junctions1 = mocker.mock()
        self.junctions2 = mocker.mock()
        self.junc1 = [400]
        self.junc2 = [700]

        juncKey = 'chr1:200'
        self.junctions1[juncKey]
        mocker.result(self.junc1)

        self.junctions2[juncKey]
        mocker.result(self.junc2)

        juncKey = 'chr1:500'
        self.junctions1[juncKey]
        mocker.throw(KeyError)

        juncKey = 'chr1:800'
        self.junctions1[juncKey]
        mocker.throw(KeyError)


        mocker.replay()
        self.sharedExons = {}
        jc.scanJunctions(self.model, self.junctions1,
                self.junctions2, self.sharedExons)

        for j in self.sharedExons:
            self.assertEqual(self.sharedExons[j], [400])
            self.assertEqual(j, 'chr1:200')

class TestBuildJunctionDict(TestCase):
    def setUp(self):
        self.mocker = Mocker()
        self.junctions = self.mocker.mock()
        self.junction1 = self.mocker.mock()
        self.junction2 = self.mocker.mock()
        self.junction3 = self.mocker.mock()


        self.junction1.start
        self.mocker.result(100)
        self.mocker.count(0, None)
        
        self.junction1.end
        self.mocker.result(200)
        self.mocker.count(0, None)

        self.junction1.chrom
        self.mocker.result('chr1')
        self.mocker.count(0, None)

        self.junction2.start
        self.mocker.result(100)
        self.mocker.count(0, None)
        
        self.junction2.end
        self.mocker.result(300)
        self.mocker.count(0, None)

        self.junction2.chrom
        self.mocker.result('chr1')
        self.mocker.count(0, None)

        self.junction3.start
        self.mocker.result(200)
        self.mocker.count(0, None)
        
        self.junction3.end
        self.mocker.result(300)
        self.mocker.count(0, None)

        self.junction3.chrom
        self.mocker.result('chr1')
        self.mocker.count(0, None)

    def test_no_overlaps(self):

        self.junctions.iteritems()
        self.mocker.generate([('chr1:100-200', self.junction1)])

        self.mocker.replay()

        self.container = jc.buildJunctionDict(self.junctions)
        self.assertEqual(self.container['chr1:100'], [200])
        self.assertEqual(len(self.container), 1)

        self.mocker.restore()
        self.mocker.verify()

    def test_overlaps(self):

        self.junctions.iteritems()
        self.mocker.generate([('chr1:100-200', self.junction1),
                            ('chr1:100-300', self.junction2)])

        self.mocker.replay()
        self.container = jc.buildJunctionDict(self.junctions)
        self.assertEqual(self.container['chr1:100'], [200, 300])
        self.assertEqual(len(self.container), 1)

        self.mocker.restore()
        self.mocker.verify()

    def test_multiple_junctions_with_overlaps(self):
        self.junctions.iteritems()
        self.mocker.generate([('chr1:100-200', self.junction1),
                            ('chr1:100-300', self.junction2),
                            ('chr1:200-300', self.junction3)],)

        self.mocker.replay()
        self.container = jc.buildJunctionDict(self.junctions)
        self.assertEqual(self.container['chr1:100'], [200, 300])
        self.assertEqual(self.container['chr1:200'], [300])
        self.assertEqual(len(self.container), 2)

        self.mocker.restore()
        self.mocker.verify()


if __name__ == '__main__':
    main()
