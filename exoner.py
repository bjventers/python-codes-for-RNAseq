#! /usr/local/bin/python
'''
    This script create exon structure by reading BLAT output to construct exon junctions.
    Requires psl-parser.py
'''

import sys
import psl_parser
import pysam

def construct(all_transcripts):
    ''' This function construct exon and exon junction from transcript model. '''
    exons = {}
    junctions = {}
    for tr in all_transcripts:
        starts = []
        ends = []
        chr = tr.attrib['tName']
        if chr not in exons:
            exons[chr] = {}
        if chr not in junctions:
            junctions[chr] = {}

        for i in range(len(tr.attrib['tStarts'])):
            start = tr.attrib['tStarts'][i] + 1 # 1-based coordinate, according to UCSC genome browser.
            size = tr.attrib['blockSizes'][i] - 1 # start position inclusive.
            end = start + size
            try:
                exons[chr][start][end] += 1
            except KeyError:
                exons[chr][start] = {end:1}
            starts.append(start)
            ends.append(end)
        for i in range(len(ends[:-1])):
            try:
                if starts[i+1] not in junctions[chr][ends[i]]:
                    junctions[chr][ends[i]].append(starts[i+1])
            except KeyError:
                junctions[chr][ends[i]] = [starts[i+1]]
    return exons, junctions

def find_SNP(chr, ex_start, ex_end, SNP_file):
    SNP = {}
    snp_file = open(SNP_file)
    for line in snp_file:
        snp_chr, pos, ref, snp, cov = line.strip().split()
        if chr == snp_chr:
            if int(pos) in range(ex_start + 10, ex_end + 10):
                return ref, geno, cov
    return None

def export_junctions(junctions):
    '''
        This function export junctions info. to standard output.
    '''
    for chr in junctions:
        for intron_start in junctions[chr]:
            intron_ends = [intron_end for intron_end in junctions[chr][intron_start]]
            intron_ends.sort()
            print '%d,' % (intron_start),
            for intron_end in intron_ends:
                print '%d,' % (intron_end),
            print

def cassette_exon(exons, junctions, samfile):
    ''' This function detects mutually exclusive exon and
        return a coverage of each exon included in each alternative splicing.
    '''
    for chr in junctions:
        for intron_start in junctions[chr]:
            if len(junctions[chr][intron_start]) > 1:
                ''' If the start site has more than one end,
                    this means there is an alternative splicing.
                '''
                intron_ends = [intron_end for intron_end in junctions[chr][intron_start]]
                intron_ends.sort()
                for intron_end in intron_ends:
                    ex_start = intron_end               # intron end site = exon start site
                    try:
                        for ex_end in exons[chr][ex_start]:
                            if len(junctions[chr][ex_end]) == 1:
                                if junctions[chr][ex_end][0] == intron_ends[-1]:
                                    coverage = len([alread for alread in samfile.fetch(chr, ex_start , ex_end)])
                                    print '%d-%d,%d,%d,%d' % (intron_start, intron_ends[-1], ex_start, ex_end, coverage),
                                    for line in open(SNP_file):
                                        snp_chr, pos, ref, gen, snp_cov = line.strip().split()
                                        if snp_chr == chr:
                                            if int(pos) in range(ex_start - 10, ex_start + 10) or int(pos) in range(ex_end - 10, ex_end + 10):
                                                print ',%s:%s:%s:%s' % (ref, gen, pos, snp_cov),
                                    print 
                    except KeyError:
                        pass

def AS_coverage(exons, junctions, samfile):
    ''' This function detects alternative splicing of each exon and
        return a coverage of each exon included in each alternative splicing.
    '''
    for chr in junctions:
        for intron_start in junctions[chr]:               # ss = intron start site
            if len(junctions[chr][intron_start]) > 1:
                ''' If the start site has more than one end,
                    this means there is an alternative splicing.
                '''
                for exon_start in exons[chr]:
                    for exon_end in exons[chr][exon_start]:
                        if exon_end == intron_start:
                            as_exons = [exon_start]
                            coverage = len([alread for alread in samfile.fetch(chr, exon_start, exon_end)])
                            print '>%d\t%d\t%d' % (exon_start, exon_end, coverage)
                            break

                for intron_end in junctions[chr][intron_start]:   # es = intron end site
                    ex_start = intron_end               # intron end site = exon start site
                    for ex_end in exons[chr][ex_start]:
                        coverage = len([alread for alread in samfile.fetch(chr, ex_start , ex_end)])
                        print '%d\t%d\t%d' % (ex_start, ex_end, coverage)
                        as_exons.append(ex_end)
                as_exons.sort()
                print '%d-%d' % (as_exons[0], as_exons[-1])

if __name__ == '__main__':

    fname = sys.argv[1]
    samfile = pysam.Samfile(sys.argv[2], 'rb')
    SNP_file = sys.argv[3]

    comment = 'track'
    all_transcripts = []
    for line in psl_parser.read(open(fname), comment):
        all_transcripts.append(line)

    exons, junctions = construct(all_transcripts)
    #AS_coverage(exons, junctions, samfile)
    #cassette_exon(exons, junctions, samfile)
    export_junctions(junctions)
