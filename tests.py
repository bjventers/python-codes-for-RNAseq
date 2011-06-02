import exon_g as exon_grouper
import psl_parser
def test_foo():
    exons = {}
    for aln in psl_parser.read(open('rap1b.psl'), 'track'):
        exon_grouper.construct(aln.attrib['tStarts'], aln.attrib['blockSizes'], aln.attrib['tName'], exons)

    multiple_junctions = [exon for exon in exons.itervalues() if len(exon.junctions) > 1]
    exon_clusters = exon_grouper.cluster(exons)
    #exon_grouper.printout_BED(exons, exon_clusters)
    print exon_clusters
    assert exon_clusters == [37086008, 37115513]

def test_cluster():
    exons = {}
    tStarts = [100, 300, 500, 700]
    blockSizes = [100, 50, 150, 200]
    exon_grouper.construct(tStarts, blockSizes, 'chr1', exons)

    multiple_junctions = [exon for exon in exons.itervalues() if len(exon.junctions) > 1]
    exon_clusters = exon_grouper.cluster(exons)
    #exon_grouper.printout_BED(exons, exon_clusters)
    print exon_clusters
    assert exon_clusters == [200]
    
