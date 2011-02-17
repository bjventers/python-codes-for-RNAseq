#! /usr/local/bin/python

'''
    This script groups all exons related to the same transcript together.
    The output is in BED format for viewing in UCSC genome browser.
    The script requires the alignment of transcript assembly from
    velvet + oases output to the referecne genome.
    The input has to be in PSL/BED format from GMAP, BLAT or this script.

    Author: Likit Preeyanon
    Email: preeyano@msu.edu
'''

import psl_parser
import sys


class Exon(object):

    def __init__(self, chr, start, end, junctions):
        self.chr = chr
        self.start = start
        self.end = end
        self.junctions = junctions
        self.cluster = []
        self.connected_exons = [(sorted(start)[0], end)]


class BED(object):

    def __init__(self, **kwargs):
        self.attrib = kwargs
        self.attrib['blockSizes'] = [int(i) for i in self.attrib['blockSizes']]
        self.attrib['tStarts'] = [int(i) for i in self.attrib['tStarts']]


def parse_bed(fp, comment):

    attrib = {}
    for line in fp:
        if line.startswith(comment):
            continue
        rows = line.strip().split()
        attrib['tName'] = rows[0]
        attrib['tStarts'] = rows[-1].split(',')[:-1]
        attrib['blockSizes'] = rows[-2].split(',')[:-1]
        if len(attrib['tStarts']) == 1:
            continue
        aln = BED(**attrib)
        yield aln

def construct(aln_obj, exons):

    for i in range(len(aln.attrib['tStarts'][:-1])):
        end = aln.attrib['tStarts'][i] + aln.attrib['blockSizes'][i]
        start = aln.attrib['tStarts'][i]
        if end in exons and aln.attrib['tName'] == exons[end].chr:
            if start not in exons[end].start:
                exons[end].start.append(start)
            junc = (aln.attrib['tStarts'][i + 1], \
                    aln.attrib['tStarts'][i + 1] + \
                    aln.attrib['blockSizes'][i + 1])

            if junc not in exons[end].junctions:
                exons[end].junctions.append(junc)
        else:
             exons[end] = Exon(aln.attrib['tName'], [start], \
                                 end, [(aln.attrib['tStarts'][i + 1], \
                                 aln.attrib['tStarts'][i + 1] + \
                                 aln.attrib['blockSizes'][i + 1])])

    last_exon_start = aln.attrib['tStarts'][-1]
    last_exon_end = last_exon_start + aln.attrib['blockSizes'][-1]

    if last_exon_end in exons and \
        aln.attrib['tName'] == exons[last_exon_end].chr:

        if last_exon_start not in exons[last_exon_end].start:
            exons[last_exon_end].start.append(last_exon_start)
    else:
        exons[last_exon_end] = Exon(aln.attrib['tName'], \
                                    [last_exon_start], last_exon_end, [])


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

    grouped.append(exon_end)

    if exons[exon_end].junctions:
        for junc in exons[exon_end].junctions:
            start, end = junc
            all_connected_exons.append((sorted(exons[end].start)[0], exons[end].end))
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
    print >> sys.stderr, 'Clustering ...'
    for num, e in enumerate(sorted(exons), start=1):
        if e in grouped:
            if num % 1000 == 0:
                print >> sys.stderr, '...', num
            continue
        else:
            '''
                All exons that are connected to this exon are saved in
                all_connected_exons.
            '''
            all_connected_exons = []  

            '''
                All connected exons will be clustered into the cluster that
                this exon belongs to.
                Note that each exon at least belong to its own cluster.
                All clusters that each exon belongs to can be looked up in
                cluster attribute of Exon object. 
            '''
            exons[e].cluster.append(e)
            add_clusters = []

            '''
                join() walks through all exons that connected together
                and add them in all_connected_exons list.
                Also, all clusters that found in cluster attribute of each
                connected exon will be added to add_clusters list.
            '''
            join(exons, exons[e].end, grouped, e, all_connected_exons, add_clusters)

            if add_clusters:
                for start, end in all_connected_exons:
                    for c in add_clusters:
                        if c not in exons[end].cluster:
                            exons[end].cluster.append(c)
                        if (start, end) not in exons[c].connected_exons:
                            exons[c].connected_exons.append((start, end))
                for c in add_clusters:
                    if c not in exons[e].cluster:
                        exons[e].cluster.append(c)

                    if (sorted(exons[e].start)[0], exons[e].end) not in exons[c].connected_exons:
                        exons[c].connected_exons.append((sorted(exons[e].start)[0], exons[e].end))

        if num % 1000 == 0:
            print >> sys.stderr, '...', num

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
            for j in connected_exons:
                if j[0] not in new_junctions:
                    new_junctions[j[0]] = j[-1]
                else:
                    if j[-1] > new_junctions[j[0]]:
                        new_junctions[j[0]] = j[-1]
        if len(new_junctions) > 1:
            chromStart = connected_exons[0][0]
            blockStarts = [j - chromStart for j in sorted(new_junctions)]
            blockSizes = [new_junctions[j] - j for j in sorted(new_junctions)]
            #blockEnds = [int(new_junctions[j]) for j in sorted(new_junctions)]

            new_blockStarts = [blockStarts[0]]  #new blockStarts
            new_blockSizes = []  #new blockSizes
            cum_size = 0  #cumulative size

            i = 0
            while i < len(blockStarts):
                current_size = blockStarts[i] + blockSizes[i]
                cum_size += blockSizes[i]
                #print '-->', blockStarts[i], current_size, cum_size
                try:
                    if blockStarts[i+1] < current_size:
                        cum_size -= blockStarts[i+1] - blockStarts[i]
                        i += 1
                    if blockStarts[i+1] == current_size:
                        i += 1
                    else:
                        new_blockSizes.append(cum_size)
                        i += 1
                        new_blockStarts.append(blockStarts[i])
                        cum_size = 0
                except IndexError:
                    if blockStarts[i] < blockStarts[i-1] + blockSizes[i-1]:
                        cum_size = (blockSizes[i] + blockSizes[i-1]) - (blockStarts[i] - blockStarts[i-1]) 
                        #print cum_size, blockSizes[i], blockSizes[i-1]
                        new_blockSizes.append(cum_size)
                    elif blockStarts[i] == blockStarts[i-1] + blockSizes[i-1]:
                        cum_size = blockSizes[i] + blockSizes[i-1]
                        new_blockSizes.append(cum_size)
                    else:
                        new_blockSizes.append(blockSizes[i])
                    break
            chromEnd = new_blockStarts[-1] + new_blockSizes[-1] + chromStart
            new_blockStarts = [str(i) for i in new_blockStarts]
            new_blockSizes = [str(i) for i in new_blockSizes]
            blockCount = len(new_blockStarts)
            print >> sys.stdout, '%s\t%d\t%d\ttest\t1000\t+\t%d\t%d\t0,0,255\t%d\t%s\t%s' % (exons[e].chr, chromStart, chromEnd, chromStart, chromEnd, blockCount, ','.join(new_blockSizes), ','.join(new_blockStarts))

if __name__ == '__main__':
    print >> sys.stderr, 'Constructing exons ...'
    exons = {}
    for aln in psl_parser.read(open(sys.argv[1]), 'track'):
        construct(aln, exons)
    print >> sys.stderr, 'total exons = %d' % len(exons)
    exon_clusters = cluster(exons)
    printout_BED(exons, exon_clusters)
