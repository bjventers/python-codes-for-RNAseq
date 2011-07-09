import sys
from pygr import seqdb
from Bio import SeqIO
from Bio.Blast.Applications import NcbiblastpCommandline
from Bio.Blast import NCBIWWW, NCBIXML

for seqRecord in SeqIO.parse('proteins.fa', 'fasta'):
    print 'doing BLAST search for', seqRecord.id
    ratios = []
    resultHandle = NCBIWWW.qblast('blastp', 'nr', seqRecord.seq)
    blastRecords = NCBIXML.parse(resultHandle)
    for blastRecord in blastRecords:
        for alignment in blastRecord.alignments:
            for hsp in alignment.hsps:
                if hsp.expect < 0.04 and 'Mus' in alignment.title: 
                    #print '----Alignment----'
                    print 'sequence:', alignment.title
                    #print 'length:', alignment.length
                    #print 'bits:', hsp.bits
                    #print 'identity:', hsp.identities
                    print 'bits/length:', hsp.bits/alignment.length
                    print 'identity/length:', hsp.identities/float(alignment.length)
                    ratios.append(hsp.bits/alignment.length)
                    '''
                    print hsp.query[0:75] + '...'
                    print hsp.match[0:75] + '...'
                    print hsp.sbjct[0:75] + '...'
                    '''
    print sorted(ratios)

#def validateORF(seqRecord):
#blastpCommand = NcbiblastpCommandline(query='proteins.fa', db='nr', out='proteinSample.xml')
