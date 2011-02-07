import pysam
import sys

'''Input file is BAM file with BAM index file'''

if len(sys.argv) > 3:
    chr = sys.argv[2]
    start = int(sys.argv[3])
    stop = int(sys.argv[4])
elif len(sys.argv) > 2:
    chr = sys.argv[2]
    start = None
    stop = None

samfile = pysam.Samfile(sys.argv[1], 'rb')

for pileupcolumn in samfile.pileup(chr, start, stop):
    print >> sys.stdout, chr, pileupcolumn.pos, pileupcolumn.pos+1, pileupcolumn.n

