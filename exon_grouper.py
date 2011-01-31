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

def join(exons, ex_end, all):
    if exons[ex_end].junctions:
        for exon in exons[ex_end].junctions:
            start, end = exon
            all.append((start, end))
            join(exons, end, all)
    else:
        return

if __name__ == '__main__':
    exons = {}
    for aln in psl_parser.read(open(sys.argv[1]), 'track'):
        construct(aln, exons)

    for e in exons:
        if len(exons[e].junctions) == 3:
            all = []
            join(exons, exons[e].end, all)
            if len(all) > len(exons[e].junctions):
                #print exons[e].chr, exons[e].start, exons[e].end, exons[e].junctions
                #print 'after grouping:\n', sorted(all)
                break
    new_all = {}
    for j in all:
        if j[0] not in new_all:
            new_all[j[0]] = j[-1]
        else:
            if j[-1] != new_all[j[0]]:
                new_all[j[0]] = j[-1]

    chromStart = sorted(new_all.keys())[0]
    blockStarts = [str(j - chromStart) for j in sorted(new_all)]
    blockSizes = [str(new_all[j] - j) for j in sorted(new_all)]
    blockEnds = [str(new_all[j] - chromStart) for j in sorted(new_all)]
    chromEnd = int([str(new_all[j]) for j in sorted(new_all)][-1])
    blockCount = len(new_all)
    print '%s\t%d\t%d\ttest\t1000\t+\t%d\t%d\t255,0,0\t%d\t%s\t%s' % (exons[e].chr, chromStart, chromEnd, chromStart, chromEnd, blockCount, ','.join(blockSizes), ','.join(blockStarts))
