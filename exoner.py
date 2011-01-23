#! /usr/local/bin/python
'''
    This script create exon structure by reading BLAT output to construct exon junctions.
    Requires psl-parser.py
'''

import sys
import psl_parser

fname = sys.argv[1]

for line in psl_parser.read(open(fname)):
    print line['qName'], line['tName'], line['qStarts'], line['tStarts']
    break
