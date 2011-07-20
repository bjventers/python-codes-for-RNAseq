'''This script blast the query file to the protein

database (nr--for internet BLAST, user specified db--for local BLAST)
and report the highest bit score:length ratio. E-value can be specified.

'''

import sys
import argparse

from Bio import SeqIO
from Bio.Blast import NCBIWWW, NCBIXML


def internetBLAST(inputFile, fileFormat='fasta', evalue=0.001):
    for seqRecord in SeqIO.parse(inputFile, fileFormat):
        print >> sys.stderr, 'Doing BLAST (internet) search for', seqRecord.id
        ratios = []
        resultHandle = NCBIWWW.qblast('blastp', 'nr', seqRecord.seq)
        blastRecords = NCBIXML.parse(resultHandle)
        for blastRecord in blastRecords:
            for alignment in blastRecord.alignments:
                for hsp in alignment.hsps:
                    if hsp.expect < evalue:
                        ratios.append(hsp.bits / alignment.length)
        if ratios:
            print >> sys.stdout, '%s\t%f' % (seqRecord.id, max(ratios))
        else:
            print >> sys.stdout, '%s\t%s' % (seqRecord.id, 'NA')


def localBLAST_gene(inputFile, evalue=0.001):
    handle = open(inputFile)
    blastRecords = NCBIXML.parse(handle)
    geneName = None
    for blastRecord in blastRecords:
        queryName = blastRecord.query.split()[0]
        chrom, newGeneID, isoformID = queryName.split('_')
        newGeneName = chrom + '_' + newGeneID
        if not geneName:
            geneName = newGeneName
            subjects = []
            ratios = []
        else:
            '''If the record is a new gene, print out

            the result of the previous gene

            '''
            if newGeneName != geneName:
                if ratios:
                    for geneName, ratio, \
                            subj in ratios:
                        print >> sys.stdout, '%s\t%.16f\t%s' % (geneName,
                                                                ratio,
                                                                subj)
                        print >> sys.stderr, '%s\t%.16f\t%s' % (geneName,
                                                                ratio,
                                                                subj)
                ratios = []
                subjects = []
                geneName = newGeneName
            else:
                '''If the record is another isoform,

                calculate the ratio and store the maximum
                ratio in the list.

                '''
                maxRatio = None
                maxRatioSubj = None
                for alignment in blastRecord.alignments:
                    for hsp in alignment.hsps:
                        subject = alignment.title.split()[0]
                        if hsp.expect < evalue:
                            ratio = hsp.bits / len(hsp.sbjct)
                            if ratio > maxRatio:
                                maxRatio = ratio
                                maxRatioSubj = subject
                if maxRatio:
                    if maxRatioSubj not in subjects:
                        '''If the subject sequence is a new sequence,

                        save the maximum ratio in the list.

                        '''
                        ratios.append((geneName, maxRatio, maxRatioSubj))
                        subjects.append(maxRatioSubj)
                    else:
                        '''If the subject sequence is existing,

                        update the ratio if the new ratio is higher.

                        '''
                        for i in range(len(ratios)):
                            q, ratio, subj = ratios[i]
                            if (subj == maxRatioSubj and maxRatio > ratio):
                                ratios[i] = (geneName, maxRatio, maxRatioSubj)

    '''Print out the result of the last record'''
    if ratios:
        for geneName, ratio, subj in ratios:
            print >> sys.stdout, '%s\t%.16f\t%s' % (geneName, ratio, subj)


def localBLAST_isoform(inputFile, evalue=0.01):
    handle = open(inputFile)
    blastRecords = NCBIXML.parse(handle)
    for blastRecord in blastRecords:
        isoform = blastRecord.query.split()[0]
        maxRatio = None
        maxRatioSubj = None
        for alignment in blastRecord.alignments:
            for hsp in alignment.hsps:
                subject = alignment.title.split()[0]
                if hsp.expect < evalue:
                    ratio = hsp.bits / len(hsp.sbjct)
                    if ratio > maxRatio:
                        maxRatio = ratio
                        maxRatioSubj = subject
        if maxRatio:
            print >> sys.stdout, '%s\t%.16f\t%s' % (isoform,
                                                    maxRatio,
                                                    maxRatioSubj)
            print >> sys.stderr, '%s\t%.16f\t%s' % (isoform,
                                                    maxRatio,
                                                    maxRatioSubj)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
            formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('--method', dest='method', default='local',
                        help='BLAST method [local, internet]')
    parser.add_argument('--evalue', dest='evalue', type=float,
                        help='e-value', default=0.001)
    parser.add_argument('--format', dest='fileFormat',
                        help='input file format', default='fasta')
    parser.add_argument('--input', dest='inputFile',
                        help='input file')
    parser.add_argument('--output', dest='outputFormat', default='gene',
                        help='output format [isoform, gene]')
    args = parser.parse_args()

    if args.method == 'local':
        if args.outputFormat == 'gene':
            localBLAST_gene(args.inputFile, args.evalue)
        else:
            localBLAST_isoform(args.inputFile, args.evalue)
    else:
        internetBLAST(args.inputFile, args.fileFormat, args.evalue)
