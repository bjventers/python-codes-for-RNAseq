#! /usr/local/bin/python

'''
    This script parse exon junctions from tophat output, junctions.bed, file
    and find alternative splicings.
    Author: Likit Preeyanon
    Email: preeyano@msu.edu
'''

import sys

try:
    for line in open(fp):

        if line.startswith('track'):
            continue

        rows = line.split()
        chrom = rows[0]
        intron_start, intron_end = rows[10].split(',')
        intron_start = int(intron_start)
        intron_end = int(intron_end)
        start = int(rows[1]) + intron_start
        end = int(rows[2]) - intron_end
        strand = rows[5]
        print chrom, start, end
        break
except IOError:
    print >> sys.stderr, 'Check the input file ...'
