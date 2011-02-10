#! /usr/local/bin/python

'''
    This script parse exon junctions from tophat output, junctions.bed,
    and find alternative splicings.
    The script handles only one reference sequence per input file.
    Use: grep -w 'reference sequence' junctions.bed > ref.junctions.bed
    to select the sequence of interest from junctions.bed.

    Author: Likit Preeyanon
    Email: preeyano@msu.edu
'''

import sys

fp1 = sys.argv[1]  #input file1: [junctions.bed]
fp2 = sys.argv[2]  #input file2: [junctions.bed]

def parse(fp):
    junctions = {}
    try:

        '''
            Parse data from BED format.
        '''

        for line in open(fp):

            if line.startswith('track'):
                continue

            rows = line.split()
            chrom = rows[0]
            name = rows[3]
            intron_start, intron_end = rows[10].split(',')
            intron_start = int(intron_start)
            intron_end = int(intron_end)
            start = int(rows[1]) + intron_start
            end = int(rows[2]) - intron_end
            strand = rows[5]

            try:
                junctions[start].append(end)
            except KeyError:
                junctions[start] = [end]
        return junctions

    except IOError:
        print >> sys.stderr, 'Check the input file ...'
        return None

if __name__=='__main__':

    junctions1 = parse(fp1)
    junctions2 = parse(fp2)

    junctions_diff= {}

    for j in junctions1:
        if j in junctions2:
            junctions_diff[j] = []
            for end in junctions2[j]:
                if end not in junctions1[j]:
                    junctions_diff[j].append(end)

    print 'junctions1 = ', len(junctions1)
    print 'junctions2 = ', len(junctions2)
    print 'diff junctions', len([j for j in junctions_diff if junctions_diff[j]])

    for n, start in enumerate(junctions_diff):
        if junctions_diff[start]:
            print start, junctions1[start], junctions_diff[start]
        if n == 50: break
    print '...'

    '''
    for start in sorted(junctions):
        if len(junctions[start]) >= 2:
            for end in junctions[start]:
                if end - start >= 10:
                    print >> sys.stdout, '%s:%d-%d' % (chrom, start, end)
    '''
