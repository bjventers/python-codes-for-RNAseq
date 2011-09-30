import unittest
from sys import stderr
import compareGeneModels as comp


class TestDb(unittest.TestCase):
    def setUp(self):
        self.a1 = comp.ExonObj('1', 100, 200)
        self.b1 = comp.ExonObj('1', 300, 400)
        self.c1 = comp.ExonObj('1', 500, 600)
        self.d1 = comp.ExonObj('1', 700, 800)
        self.e1 = comp.ExonObj('1', 900, 1000)
        self.f1 = comp.ExonObj('1', 1100, 1200)
        self.g1 = comp.ExonObj('1', 1300, 1400)

        self.a2 = comp.ExonObj('1', 100, 200)
        self.b2 = comp.ExonObj('1', 300, 400)
        self.c2 = comp.ExonObj('1', 500, 600)
        self.d2 = comp.ExonObj('1', 700, 800)
        self.e2 = comp.ExonObj('1', 900, 1000)
        self.f2 = comp.ExonObj('1', 1100, 1200)
        self.g2 = comp.ExonObj('1', 1300, 1400)

        self.a3 = comp.ExonObj('1', 100, 200)
        self.b3 = comp.ExonObj('1', 300, 400)
        self.c3 = comp.ExonObj('1', 500, 600)
        self.d3 = comp.ExonObj('1', 700, 800)
        self.e3 = comp.ExonObj('1', 900, 1000)
        self.f3 = comp.ExonObj('1', 1100, 1200)
        self.g3 = comp.ExonObj('1', 1300, 1400)

    def testAddExonOneCluster(self):
        db = {}
        clusters = []
        exons1 = [self.a1, self.b1, self.c1, self.d1]
        exons2 = [self.a2, self.b2, self.c2, self.e2]

        comp.addExon(db, clusters, exons1)
        self.assertEqual(len(db), 4)

        comp.addExon(db, clusters, exons2)
        self.assertEqual(len(db), 5)

        self.assertEqual(len(clusters), 1)

        cl = clusters[0]

        self.assertEqual(len(cl.neighbors(self.c1.coord)), 2)
        self.assertIn(self.d1.coord, cl.nodes())
        self.assertIn(self.e1.coord, cl.nodes())
        self.assertEqual(len(cl.predecessors(self.c1.coord)), 1)
        self.assertIn(self.b1.coord, cl.predecessors(self.c1.coord))
        self.assertIn(self.a1.coord, cl.predecessors(self.b1.coord))
        self.assertEqual(cl.predecessors(self.a1.coord), [])
        self.assertEqual(cl.successors(self.d1.coord), [])

    def testAddExonTwoClusters(self):
        db = {}
        clusters = []
        exons1 = [self.a1, self.b1, self.c1, self.d1]
        exons2 = [self.a2, self.b2, self.c2, self.e2]
        exons3 = [self.g3, self.f3]

        comp.addExon(db, clusters, exons1)
        self.assertEqual(len(db), 4)

        comp.addExon(db, clusters, exons2)
        self.assertEqual(len(db), 5)

        comp.addExon(db, clusters, exons3)
        self.assertEqual(len(db), 7)

        self.assertEqual(len(clusters), 2)

        self.assertEqual(len(clusters[0].nodes()), 5)
        self.assertEqual(len(clusters[1].nodes()), 2)


class TestWalkDown(unittest.TestCase):
    def setUp(self):
        self.a1 = comp.ExonObj('1', 100, 200)
        self.b1 = comp.ExonObj('1', 300, 400)
        self.c1 = comp.ExonObj('1', 500, 600)
        self.d1 = comp.ExonObj('1', 700, 800)
        self.e1 = comp.ExonObj('1', 900, 1000)
        self.f1 = comp.ExonObj('1', 1100, 1200)
        self.g1 = comp.ExonObj('1', 1300, 1400)

        self.a2 = comp.ExonObj('1', 100, 200)
        self.b2 = comp.ExonObj('1', 300, 400)
        self.c2 = comp.ExonObj('1', 500, 600)
        self.d2 = comp.ExonObj('1', 700, 800)
        self.e2 = comp.ExonObj('1', 900, 1000)
        self.f2 = comp.ExonObj('1', 1100, 1200)
        self.g2 = comp.ExonObj('1', 1300, 1400)

        self.a3 = comp.ExonObj('1', 100, 200)
        self.b3 = comp.ExonObj('1', 300, 400)
        self.c3 = comp.ExonObj('1', 500, 600)
        self.d3 = comp.ExonObj('1', 700, 800)
        self.e3 = comp.ExonObj('1', 900, 1000)
        self.f3 = comp.ExonObj('1', 1100, 1200)
        self.g3 = comp.ExonObj('1', 1300, 1400)

    def testSinglePath(self):
        db = {}
        clusters = []
        exons1 = [self.a1, self.b1, self.c1, self.d1, self.e1]
        comp.addExon(db, clusters, exons1)

        path = []
        allPath = []

        comp.walkDown(self.a1.coord, path, allPath, self.a1.cluster) 
        self.assertEqual(len(allPath), 1)
        self.assertEqual(allPath, [[str(self.b1), str(self.c1),
                                    str(self.d1), str(self.e1)]])

    def testTwoPaths(self):
        db = {}
        clusters = []
        exons1 = [self.a1, self.b1, self.c1, self.d1]
        exons2 = [self.a2, self.b2, self.c2, self.e2]
        comp.addExon(db, clusters, exons1)
        comp.addExon(db, clusters, exons2)

        path = []
        allPath = []

        comp.walkDown(self.c2.coord, path, allPath, self.c2.cluster)
        self.assertEqual(len(allPath), 2)

        path = []
        allPath = []

        comp.walkDown(self.a2.coord, path, allPath, self.a2.cluster)
        self.assertEqual(len(allPath), 2)
        self.assertEqual(allPath[0], [str(self.b1), str(self.c1), str(self.d1)])
        self.assertEqual(allPath[1], [str(self.b2), str(self.c2), str(self.e2)])

    def testTwoPathsTwoBubbles(self):
        db = {}
        clusters = []
        exons1 = [self.a1, self.b1, self.d1, self.f1, self.g1]
        exons2 = [self.a2, self.c2, self.d2, self.e2, self.g2]
        comp.addExon(db, clusters, exons1)
        comp.addExon(db, clusters, exons2)

        path = []
        allPath = []

        comp.walkDown(self.a1.coord, path, allPath, self.a1.cluster)
        self.assertEqual(len(allPath), 4)

        self.assertEqual(allPath[0], [str(self.b1), str(self.d1),
                                        str(self.f1), str(self.g1)])
        self.assertEqual(allPath[1], [str(self.b1), str(self.d1),
                                        str(self.e2), str(self.g1)])
        self.assertEqual(allPath[2], [str(self.c2), str(self.d1),
                                        str(self.f1), str(self.g1)])
        self.assertEqual(allPath[3], [str(self.c2), str(self.d1),
                                        str(self.e2), str(self.g1)])


class TestgetPath(unittest.TestCase):
    def setUp(self):
        self.a1 = comp.ExonObj('1', 100, 200)
        self.b1 = comp.ExonObj('1', 300, 400)
        self.c1 = comp.ExonObj('1', 500, 600)
        self.d1 = comp.ExonObj('1', 700, 800)
        self.e1 = comp.ExonObj('1', 900, 1000)
        self.f1 = comp.ExonObj('1', 1100, 1200)
        self.g1 = comp.ExonObj('1', 1300, 1400)

        self.a2 = comp.ExonObj('1', 100, 200)
        self.b2 = comp.ExonObj('1', 300, 400)
        self.c2 = comp.ExonObj('1', 500, 600)
        self.d2 = comp.ExonObj('1', 700, 800)
        self.e2 = comp.ExonObj('1', 900, 1000)
        self.f2 = comp.ExonObj('1', 1100, 1200)
        self.g2 = comp.ExonObj('1', 1300, 1400)

        self.a3 = comp.ExonObj('1', 100, 200)
        self.b3 = comp.ExonObj('1', 300, 400)
        self.c3 = comp.ExonObj('1', 500, 600)
        self.d3 = comp.ExonObj('1', 700, 800)
        self.e3 = comp.ExonObj('1', 900, 1000)
        self.f3 = comp.ExonObj('1', 1100, 1200)
        self.g3 = comp.ExonObj('1', 1300, 1400)

    def testSinglePath(self):
        db = {}
        clusters = []
        exons1 = [self.a1, self.b1, self.c1, self.d1, self.e1, self.f1, self.g1]
        comp.addExon(db, clusters, exons1)

        allPath = comp.getPath(self.e1.cluster)
        self.assertEqual(len(allPath), 1)
        self.assertEqual(len(allPath[0]), 7)

    def testTwoPaths(self):
        db = {}
        clusters = []
        exons1 = [self.a1, self.b1, self.c1, self.d1, self.f1]
        exons2 = [self.a2, self.b2, self.c2, self.e2, self.g2]

        comp.addExon(db, clusters, exons1)
        comp.addExon(db, clusters, exons2)

        allPath = comp.getPath(self.a1.cluster)

        self.assertEqual(len(allPath), 2)
        self.assertEqual(len(allPath[0]), 5)
        self.assertEqual(len(allPath[1]), 5)

        allPath = comp.getPath(self.c1.cluster)

        self.assertEqual(len(allPath), 2)
        self.assertEqual(len(allPath[0]), 5)
        self.assertEqual(len(allPath[1]), 5)

        allPath = comp.getPath(self.d1.cluster)

        self.assertEqual(len(allPath), 2)
        self.assertEqual(len(allPath[0]), 5)

    def testTwoPathsTwoBubbles(self):
        db = {}
        clusters = []
        exons1 = [self.a1, self.b1, self.d1, self.f1, self.g1]
        exons2 = [self.a2, self.c2, self.d2, self.e2, self.g2]
        comp.addExon(db, clusters, exons1)
        comp.addExon(db, clusters, exons2)

        allPath = comp.getPath(self.d1.cluster) 
        self.assertEqual(len(allPath), 4)

    def testTwoPathsTwoRoots(self):
        db = {}
        clusters = []
        exons1 = [self.a1, self.c1, self.d1, self.f1]
        exons2 = [self.b2, self.c2, self.e2, self.g2]

        comp.addExon(db, clusters, exons1)
        comp.addExon(db, clusters, exons2)

        allPath = comp.getPath(self.a1.cluster)

        self.assertEqual(len(allPath), 4)
        self.assertEqual(len(allPath[0]), 4)
        self.assertEqual(len(allPath[1]), 4)

        allPath = comp.getPath(self.b2.cluster)

        self.assertEqual(len(allPath), 4)
        self.assertEqual(len(allPath[0]), 4)
        self.assertEqual(len(allPath[1]), 4)

        allPath = comp.getPath(self.c1.cluster)

        self.assertEqual(len(allPath), 4)
        self.assertEqual(len(allPath[0]), 4)
        self.assertEqual(len(allPath[1]), 4)

        allPath = comp.getPath(self.d1.cluster)

        self.assertEqual(len(allPath), 4)
        self.assertEqual(len(allPath[0]), 4)
        self.assertEqual(len(allPath[1]), 4)


class TestFindPathDiff(unittest.TestCase):
    def setUp(self):
        self.a1 = comp.ExonObj('1', 100, 200)
        self.b1 = comp.ExonObj('1', 300, 400)
        self.c1 = comp.ExonObj('1', 500, 600)
        self.d1 = comp.ExonObj('1', 700, 800)
        self.e1 = comp.ExonObj('1', 900, 1000)
        self.f1 = comp.ExonObj('1', 1100, 1200)
        self.g1 = comp.ExonObj('1', 1300, 1400)

        self.a2 = comp.ExonObj('1', 100, 200)
        self.b2 = comp.ExonObj('1', 300, 400)
        self.c2 = comp.ExonObj('1', 500, 600)
        self.d2 = comp.ExonObj('1', 700, 800)
        self.e2 = comp.ExonObj('1', 900, 1000)
        self.f2 = comp.ExonObj('1', 1100, 1200)
        self.g2 = comp.ExonObj('1', 1300, 1400)

        self.a3 = comp.ExonObj('1', 100, 200)
        self.b3 = comp.ExonObj('1', 300, 400)
        self.c3 = comp.ExonObj('1', 500, 600)
        self.d3 = comp.ExonObj('1', 700, 800)
        self.e3 = comp.ExonObj('1', 900, 1000)
        self.f3 = comp.ExonObj('1', 1100, 1200)
        self.g3 = comp.ExonObj('1', 1300, 1400)

    def TestOnePathDiff(self):
        db1 = {}
        db2 = {}
        clusters = []

        exons1 = [self.a1, self.c1, self.d1, self.e1, self.f1]
        exons2 = [self.a2, self.c2, self.d2, self.e2, self.f2, self.g2]

        comp.addExon(db1, clusters, exons1)
        comp.addExon(db2, clusters, exons2)

        self.assertEqual(len(db1), 5)
        self.assertEqual(len(db2), 6)
        
        self.assertEqual(len(db1[self.a1.coord].cluster.nodes()), 5)
        self.assertEqual(len(db2[self.a1.coord].cluster.nodes()), 6)

        missingPaths, alteredPaths = comp.findPathDiff(db2, db1)

        self.assertEqual(len(missingPaths), 0)
        self.assertEqual(len(alteredPaths), 1)

    def TestTwoPathDiff(self):
        db1 = {}
        db2 = {}
        clusters = []

        exons1 = [self.a1, self.b1, self.c1, self.d1, self.f1]
        exons2 = [self.a2, self.b2, self.c2, self.d2, self.e2, self.f2]
        exons3 = [self.a3, self.b3, self.c3, self.d3, self.e3, self.g3]

        comp.addExon(db1, clusters, exons1)
        comp.addExon(db2, clusters, exons2)
        comp.addExon(db2, clusters, exons3)

        missingPaths, alteredPaths = comp.findPathDiff(db2, db1)

        self.assertEqual(len(missingPaths), 0)
        self.assertEqual(len(alteredPaths), 2)

    def TestOneMissingPath(self):
        db1 = {}
        db2 = {}
        clusters = []

        exons1 = [self.g1]
        exons2 = [self.a2, self.b2, self.c2, self.d2, self.e2, self.f2]

        comp.addExon(db1, clusters, exons1)
        comp.addExon(db2, clusters, exons2)

        missingPaths, alteredPaths = comp.findPathDiff(db2, db1)

        self.assertEqual(len(missingPaths), 1)
        self.assertEqual(len(alteredPaths), 0)


class TestFindExtendedEnds(unittest.TestCase):
    def setUp(self):
        self.aa1 = comp.ExonObj('1', 50, 200)
        self.aa2 = comp.ExonObj('1', 100, 200)
        self.b1 = comp.ExonObj('1', 300, 400)
        self.c1 = comp.ExonObj('1', 500, 600)
        self.d1 = comp.ExonObj('1', 700, 800)
        self.e1 = comp.ExonObj('1', 900, 1000)
        self.f1 = comp.ExonObj('1', 1100, 1200)
        self.gg1 = comp.ExonObj('1', 1300, 1400)
        self.gg2 = comp.ExonObj('1', 1300, 1450)

        self.a2 = comp.ExonObj('1', 100, 200)
        self.b2 = comp.ExonObj('1', 300, 400)
        self.c2 = comp.ExonObj('1', 500, 600)
        self.d2 = comp.ExonObj('1', 700, 800)
        self.e2 = comp.ExonObj('1', 900, 1000)
        self.f2 = comp.ExonObj('1', 1100, 1200)
        self.g2 = comp.ExonObj('1', 1300, 1400)

        self.a3 = comp.ExonObj('1', 100, 200)
        self.b3 = comp.ExonObj('1', 300, 400)
        self.c3 = comp.ExonObj('1', 500, 600)
        self.d3 = comp.ExonObj('1', 700, 800)
        self.e3 = comp.ExonObj('1', 900, 1000)
        self.f3 = comp.ExonObj('1', 1100, 1200)
        self.g3 = comp.ExonObj('1', 1300, 1400)

    def TestLeftExtension(self):
        db1 = {}
        db2 = {}
        clusters = []

        exons1 = [self.aa1, self.c1, self.d1, self.e1, self.f1]
        exons2 = [self.aa2, self.c2, self.d2, self.e2, self.f2, self.gg2]

        comp.addExon(db1, clusters, exons1)
        comp.addExon(db2, clusters, exons2)

        missingPaths, alteredPaths = comp.findPathDiff(db2, db1)

        self.assertEqual(len(missingPaths), 0)
        self.assertEqual(len(alteredPaths), 2)

        cluster1 = self.aa1.cluster
        cluster2 = self.aa2.cluster

        leftExtension, rightExtension = comp.findExtendedEnd(cluster1,
                                                        cluster2,
                                                            db1,
                                                            db2)
        self.assertEqual(len(leftExtension), 1)
        self.assertEqual(len(rightExtension), 0)
