'''This script makes a custom track in BED format from a list of

exons. 

Usage makeBED.py [chromosome name] [strand] [chromosome start] \
                [exon1 start,exon1 end] [...]

'''

import sys

blockStarts = []
blockSizes = []
chromName = sys.argv[1]
strand = sys.argv[2]
chromStart = int(sys.argv[3])
name = 'model_transcript'

print(sys.argv[4:])
for exon in sys.argv[4:]:
    start, end = [int(x) for x in exon.split(',')]
    blockStarts.append(start - chromStart)
    blockSizes.append(end - start)

chromEnd = end

print('%s\t%d\t%d\t%s\t%d\t%s\t%d\t%d\t%s\t%d\t%s\t%s' %
        (chromName,
        chromStart,
        chromEnd,
        name,
        1000,
        strand,
        chromStart,
        chromEnd,
        '255,0,0',
        len(blockStarts),
        ','.join([str(s) for s in blockSizes]),
        ','.join([str(s) for s in blockStarts])))
