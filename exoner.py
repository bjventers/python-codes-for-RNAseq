#! /usr/local/bin/python
'''
    This script create exon structure by reading BLAT output to construct exon junctions.
    Requires psl-parser.py
'''

import sys
import psl_parser

fname = sys.argv[1]
comment = 'track'
for line in psl_parser.read(open(fname), comment):
    print line.attrib['qName'], line.attrib['tName'], line.attrib['qStarts'], line.attrib['tStarts']
    break
