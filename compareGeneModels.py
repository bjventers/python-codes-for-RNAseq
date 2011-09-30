import sys
import csv
import networkx as nx
import matplotlib.pyplot as plt

from sys import stderr, stdout

class ExonObj(object):
    def __init__(self, chrom, start, end):
        self.chrom = chrom
        self.start = start
        self.end = end
        self.cluster = None
        self.coord = self.__str__()

    def __str__(self):
        return '%s:%d-%d' % (self.chrom, self.start, self.end)


def parseBed(filename):
    '''Reads BED file and returns a group of exons
    of a transcript at a time.

    '''

    with open(filename) as fp:
        for row in csv.reader(fp, dialect='excel-tab'):
            exons = []
            chrom = row[0]
            chromStart = int(row[1])

            '''Get all exons except terminal ones.'''
            blockStarts = [int(i) for i in row[-1].split(',')]
            blockSizes = [int(i) for i in row[-2].split(',')]

            if len(blockStarts) == 1:
                continue

            for i in range(len(blockStarts)):
                start = chromStart + blockStarts[i]
                end = start + blockSizes[i]
                exons.append(ExonObj(chrom, start, end))

            yield exons


def addExon(db, clusters, exons):
    '''
        Adds exons to a given dictionary.
        db : a dictionary to which exons will be added.
        clusters : a list to which a graph of connected
        exons will be added.

    '''

    '''Create a new cluster for all exons.'''
    newCluster = nx.DiGraph()
    newCluster.add_path([e.coord for e in exons])

    '''Check if any exon belongs to an existing cluster.'''
    exist = False
    for exonObj in exons:
        try:
            existCluster = db[exonObj.coord].cluster
        except KeyError:
            pass
        else:
            exist = True
            break

    '''If no existing cluster is found,
    create a new one with all new exons.

    '''
    if not exist:
        clusters.append(newCluster)
        for exonObj in exons:
            exonObj.cluster = newCluster
            db[exonObj.coord] = exonObj 
    else:
        ''' If an existing cluster is found,
        add all exons to it.

        '''

        existCluster.add_edges_from(newCluster.edges())

        for exonObj in exons:
            exonObj.cluster = existCluster
            db[exonObj.coord] = exonObj


def walkDown(exonCoord, path, allPath, cluster):
    '''Returns all downstream exons from a given exon.
    
    '''

    if cluster.successors(exonCoord) == []:
        allPath.append(path[:])
        return
    else:
        for nex in cluster.successors(exonCoord):
            if nex not in path:
                path.append(nex)

            walkDown(nex, path, allPath, cluster)

            path.pop()


def getPath(cluster):
    '''Returns all paths of a given cluster.'''

    roots = [node for node in cluster.nodes() \
                    if not cluster.predecessors(node)]
    allPaths = []

    for root in roots:
        path = [root]
        walkDown(root, path, allPaths, cluster)

    return allPaths


def findExtendedEnd(cluster1, cluster2, db1, db2):

    '''Returns extended left and right ends of the
    transcript.

    '''

    leftExtension = []
    rightExtension = []

    leftNodes1 = [n for n in cluster1.nodes() if not cluster1.predecessors(n)]
    leftNodes2 = [n for n in cluster2.nodes() if not cluster2.predecessors(n)]

    for n1 in leftNodes1:
        for n2 in leftNodes2:
            if set(cluster1.neighbors(n1)). \
                    intersection(set(cluster2.neighbors(n2))):
                exon1 = db1[n1]
                exon2 = db2[n2]

                if exon1.end == exon2.end:
                    if exon1.start > exon2.start:
                        leftExtension.append(n1)

    rightNodes1 = [n for n in cluster1.nodes() if not cluster1.successors(n)]
    rightNodes2 = [n for n in cluster2.nodes() if not cluster2.successors(n)]

    for n1 in rightNodes1:
        for n2 in rightNodes2:
            if set(cluster1.predecessors(n1)). \
                    intersection(set(cluster2.predecessors(n2))):
                exon1 = db1[n1]
                exon2 = db2[n2]

                if exon1.start == exon2.start:
                    if exon1.end < exon2.end:
                        rightExtension.append(n1)

    return leftExtension, rightExtension

def findPathDiff(db1, db2):
    '''Returns paths in db1 that are not in db2.'''

    checked = []
    missingPaths = []
    alteredPaths = []

    for exon in db1:
        if exon in checked:
            continue
        else:
            cl1 = db1[exon].cluster
            try:
                cl2 = db2[exon].cluster
            except KeyError:
                pass
            else:
                '''Find edges that in cluster1 not in cluster2.
                Then add paths that contains those edges to a
                list of altered paths.

                '''
                newCluster = nx.DiGraph()
                newCluster.add_edges_from([e for e in cl1.edges()\
                                            if e not in cl2.edges()])
                for node in newCluster.nodes():
                    for path in getPath(cl1):
                        if node in path and path not in alteredPaths:
                            alteredPaths.append(path)

                checked += cl1.nodes()

    for exon in db1:
        if exon not in checked:
            '''Add all paths that contains the exon
            to a list of missing paths.
            
            '''
            paths = getPath(cl1)
            missingPaths += [p for p in paths if exon in p]
            checked += cl1.nodes()

    return missingPaths, alteredPaths


def main(file1, file2):
    db1 = {}  # contains pairs of exons and their clusters.
    db2 = {}
    clusters1 = []  # contains a list of exon clusters.
    clusters2 = []
    
    for exons in parseBed(file1):
        addExon(db1, clusters1, exons)  # add exons to a database

    for exons in parseBed(file2):
        addExon(db2, clusters2, exons)

    missingPaths, alteredPaths = findPathDiff(db1, db2)

    print 'Missing paths total %d' % len(missingPaths)
    print 'Different paths total %d ' % len(alteredPaths)

    for mp in missingPaths:
        print mp

    print '\n\n'
    for ap in alteredPaths:
        print ap


if __name__=='__main__':
    main(sys.argv[1], sys.argv[2])
