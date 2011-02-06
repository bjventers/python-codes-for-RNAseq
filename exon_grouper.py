#! /usr/local/bin/python
'''
    This script groups all exons related to the same transcript together.
    The output is in PSL format for viewing in UCSC genome browser.
    The script requires the alignment of transcript assembly from oases program to the referecne genome.
    The input has to be in PSL format from GMAP or BLAT program.
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

def construct(aln_obj, exons):
    for i in range(len(aln.attrib['tStarts'][:-1])):
        end = aln.attrib['tStarts'][i] + aln.attrib['blockSizes'][i]
        start = aln.attrib['tStarts'][i]
        if end in exons and aln.attrib['tName'] == exons[end].chr:
            if start not in exons[end].start:
                exons[end].start.append(start)
            junc = (aln.attrib['tStarts'][i+1], aln.attrib['tStarts'][i+1] + aln.attrib['blockSizes'][i+1])
            if junc not in exons[end].junctions:
                exons[end].junctions.append(junc)
        else:
             exons[end] = Exon(aln.attrib['tName'], [start], end, [(aln.attrib['tStarts'][i+1], aln.attrib['tStarts'][i+1] + aln.attrib['blockSizes'][i+1])])
    last_exon_start = aln.attrib['tStarts'][-1]
    last_exon_end = last_exon_start + aln.attrib['blockSizes'][-1]
    if last_exon_end in exons and aln.attrib['tName'] == exons[last_exon_end].chr:
        if last_exon_start not in exons[last_exon_end].start: exons[last_exon_end].start.append(last_exon_start)
    else:
        exons[last_exon_end] = Exon(aln.attrib['tName'], [last_exon_start], last_exon_end, [])

def join(exons, exon_end, grouped, new_cluster, all_connected_exons, add_clusters):
    grouped.append(exon_end) # grouped exons will not be used to from a new cluster
    if exons[exon_end].junctions:
        for junc in exons[exon_end].junctions:
            start, end = junc
            all_connected_exons.append((sorted(exons[end].start)[0], exons[end].end))
            if new_cluster not in exons[end].cluster:
                exons[end].cluster.append(new_cluster) # add new cluster to the exon's cluster list.
            for c in exons[end].cluster:
                if c not in add_clusters:
                    add_clusters.append(c)
            join(exons, end, grouped, new_cluster, all_connected_exons, add_clusters)
    else:
        return

def cluster(exons):
    grouped = []
    exon_clusters = []
    print >> sys.stderr, 'Clustering ...'
    for e in sorted(exons):
        if e in grouped: 
            continue
        else:
            all_connected_exons = []
            add_clusters = []
            exons[e].cluster.append(e) # the exon belongs to its cluster
            join(exons, exons[e].end, grouped, e, all_connected_exons, add_clusters)
            if add_clusters:
                for start,end in all_connected_exons:
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
            if len(exons[e].cluster) == 1: 
                '''
                    If the cluster does not contain any exon that shared with other clusters,
                    add it to a group of unique cluster.
                '''
                exon_clusters.append(e)
    return exon_clusters

def printout_PSL(exons, exon_clusters):
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
            blockStarts = [str(j - chromStart) for j in sorted(new_junctions)]
            #blockStarts.insert(0,'0')
            blockSizes = [str(new_junctions[j] - j) for j in sorted(new_junctions)]
            #blockSizes.insert(0,str(exons[e].end-chromStart))
            blockEnds = [int(new_junctions[j]) for j in sorted(new_junctions)]
            chromEnd = blockEnds[-1]
            blockCount = len(new_junctions)
            print >> sys.stdout, '%s\t%d\t%d\ttest\t1000\t+\t%d\t%d\t0,0,255\t%d\t%s\t%s' % (exons[e].chr, chromStart, chromEnd, chromStart, chromEnd, blockCount, ','.join(blockSizes), ','.join(blockStarts))

if __name__ == '__main__':
    print >> sys.stderr, 'Constructing exons ...'
    exons = {}
    for aln in psl_parser.read(open(sys.argv[1]), 'track'):
        construct(aln, exons)
    print >> sys.stderr, 'total exons = %d' % len(exons)
    exon_clusters = cluster(exons)
    printout_PSL(exons, exon_clusters)
