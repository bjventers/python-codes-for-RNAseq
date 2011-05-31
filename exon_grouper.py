#! /usr/local/bin/python

'''
    This script groups all exons from the same transcript together.
    The output is in BED format which can be visualized in UCSC genome browser.
    The script requires the alignment of transcript assembly from
    velvet + oases output to the referecne genome.
    The input has to be in PSL format from GMAP, BLAT or this script.

    Usage: python exon_grouper.py [transcripts.psl]
    The output is written in a standard output.

    Author: Likit Preeyanon
    Email: preeyano@msu.edu
'''

import psl_parser
import sys
import time

class Exon(object):

    def __init__(self, ref, start, end, junctions):
        self.ref = ref 
        self.start = start
        self.end = end
        self.junctions = junctions
        self.cluster = []
        self.connected_exons = [(start, end)]

def construct(aln_obj, exons):
    for i in range(len(aln.attrib['tStarts'][:-1])):
        end = aln.attrib['tStarts'][i] + aln.attrib['blockSizes'][i]
        start = aln.attrib['tStarts'][i]

        if end in exons and aln.attrib['tName'] == exons[end].ref:

            if start < exons[end].start: exons[end].start = start

            junc_start = aln.attrib['tStarts'][i+1]
            junc_end = aln.attrib['tStarts'][i+1] + aln.attrib['blockSizes'][i+1]

            if (junc_start, junc_end) not in exons[end].junctions:
                exons[end].junctions.append((junc_start, junc_end))

        else:
            junc_start = aln.attrib['tStarts'][i+1]
            junc_end = aln.attrib['tStarts'][i+1] + aln.attrib['blockSizes'][i+1]
            junc = [(junc_start, junc_end)]
            exons[end] = Exon(aln.attrib['tName'], start, end, junc)

    last_exon_start = aln.attrib['tStarts'][-1]
    last_exon_end = last_exon_start + aln.attrib['blockSizes'][-1]

    if last_exon_end in exons and aln.attrib['tName'] == \
        exons[last_exon_end].ref:
        if last_exon_start < exons[last_exon_end].start:
            exons[last_exon_end].start = last_exon_start
    else:
        exons[last_exon_end] = Exon(aln.attrib['tName'], \
                                    last_exon_start, last_exon_end, [])

def join(exons, exon_end, grouped, new_cluster, \
            all_connected_exons, add_clusters):
    '''
        This function walks through all exons connected to the starting
        exon.
        Parameters:
            exons: dictionary containing Exon object.
            exon_end: key(end) of the starting exon.
            grouped: a list containing exons that already clustered.
            new_cluster: a key(end) of a starting exon to be used as a
            cluster name.
            all_connected_exons: a list to be added all exons connected to
            the starting exon.
            add_clusters: a list to be added all clusters that each exon
            belongs to. 
    '''

    if exon_end not in grouped:
        grouped.append(exon_end)

    if exons[exon_end].junctions:
        for start, end in exons[exon_end].junctions:
            all_connected_exons.append((exons[end].start, exons[end].end))
            if new_cluster not in exons[end].cluster:
                exons[end].cluster.append(new_cluster)  # add new cluster to the exon's cluster list.
            for c in exons[end].cluster:
                if c not in add_clusters:
                    add_clusters.append(c)
            join(exons, end, grouped, new_cluster, all_connected_exons, add_clusters)
    else:
        return

def cluster(exons):

    grouped = []
    exon_clusters = []
    n = 0
    now = time.time()
    print >> sys.stderr, 'Clustering ...'
    for num, e in enumerate(sorted(exons)):
        if e in grouped:
            if num % 1000 == 0:
                print >> sys.stderr, '...', num,
                print >> sys.stderr, time.time() - now, len(grouped)
                now = time.time()
            continue
        else:
            '''
                All exons that are connected to this exon are stored in
                all_connected_exons.
            '''
            all_connected_exons = []  

            '''
                All connected exons will be put into the cluster that
                this exon belongs to.
                Note that each exon belongs to its own cluster.
                All clusters that each exon belongs to can be looked up in
                cluster attribute of Exon object. 
            '''
            exons[e].cluster.append(e)
            new_clusters = []

            '''
                join() walks through all exons that connected together
                and add them in all_connected_exons list.
                Also, all clusters that found in cluster attribute of each
                connected exon will be added to new_clusters list.
            '''
            join(exons, exons[e].end, grouped, e, all_connected_exons, new_clusters)
            #print >> sys.stderr, exons[e].start, exons[e].end, all_connected_exons

            if new_clusters:
                for start, end in all_connected_exons:
                    for nc in new_clusters:
                        if nc not in exons[end].cluster:
                            exons[end].cluster.append(nc)
                        if (start, end) not in exons[nc].connected_exons:
                            exons[nc].connected_exons.append((start, end))
                for nc in new_clusters:
                    if nc not in exons[e].cluster:
                        exons[e].cluster.append(nc)

                    if (exons[e].start, exons[e].end) not in exons[nc].connected_exons:
                        exons[nc].connected_exons.append((exons[e].start, exons[e].end))

        if num % 1000 == 0:
            print >> sys.stderr, '...', num, time.time() - now, len(grouped)
            now = time.time()

    for e in sorted(exons):
        if len(exons[e].cluster) == 1 and len(exons[e].connected_exons) >= 2:
            exon_clusters.append(e)
    print >> sys.stderr, 'total clusters =', len(exon_clusters)
    return exon_clusters

def printout_BED(exons, exon_clusters):

    for e in exon_clusters:
        new_junctions = {}
        connected_exons = sorted(exons[e].connected_exons)
        if connected_exons:
            for start, end in connected_exons:
                if start not in new_junctions:
                    new_junctions[start] = end
                else:
                    if end > new_junctions[start]:
                        new_junctions[start] = end

        if new_junctions:
            chromStart = connected_exons[0][0]
            blockStarts = [j - chromStart for j in sorted(new_junctions)]
            blockSizes = [new_junctions[j] - j for j in sorted(new_junctions)]
            #blockEnds = [int(new_junctions[j]) for j in sorted(new_junctions)]

            chromEnd = blockStarts[-1] + blockSizes[-1] + chromStart
            blockCount = len(blockStarts)
            new_blockStarts = [str(i) for i in blockStarts]
            new_blockSizes = [str(i) for i in blockSizes]

            print >> sys.stdout,'%s\t%d\t%d\ttest\t1000\t+\t%d\t%d\t0, \
                                    0,255\t%d\t%s\t%s' % (exons[e].ref, 
                                                            chromStart,
                                                            chromEnd, 
                                                            chromStart, 
                                                            chromEnd,
                                                            blockCount,
                                                            ','.join(new_blockSizes),
                                                            ','.join(new_blockStarts))

if __name__ == '__main__':
    print >> sys.stderr, 'Constructing exons ...'
    exons = {}
    for aln in psl_parser.read(open(sys.argv[1]), 'track'):
        construct(aln, exons)
    print >> sys.stderr, 'total exons = %d' % len(exons)
    multiple_junctions = [exon for exon in exons.itervalues() if len(exon.junctions) > 1]
    print >> sys.stderr, 'total multiple junctions', len(multiple_junctions)
    exon_clusters = cluster(exons)
    printout_BED(exons, exon_clusters)
