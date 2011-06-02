import exon_grouper as eg
import sys
def testClusters1():
    exons = {}
    tStartsSets = [[100, 300, 500, 700],
                    [150, 300,500, 700,800],
                    ]
    blockSizesSets = [[100, 50, 120, 200],
                        [50, 50, 150, 50, 100],
                        ]
    tName = 'chr1'
    for i in range(len(tStartsSets)):
        eg.construct(tName, tStartsSets[i], blockSizesSets[i], exons)

    exonClusters, clusterReferences = eg.cluster(exons)
    geneModels = eg.buildGeneModels(exons, exonClusters, clusterReferences)
    starts, sizes = eg.printBed(geneModels)
    print >> sys.stderr, 'sizes = ', sizes
    print >> sys.stderr, 'starts = ', starts
    assert sizes == ['100', '50', '120', '50', '100']
    assert starts == ['0', '200', '400','600', '700']

def testClusters2():
    exons = {}
    tStartsSets = [[100, 300, 500, 700],
                    #[150, 300,500, 700,800],
                    [50, 150, 300, 700],
                    #[50, 300, 500],
                    #[45, 100, 300, 500, 700, 800],
                    ]
    blockSizesSets = [[100, 50, 120, 200],
                        #[50, 50, 150, 50, 100],
                        [10, 50, 320, 200],
                        #[150, 50, 50],
                        #[15, 100, 50, 120, 50, 100],
                        ]
    tName = 'chr1'
    for i in range(len(tStartsSets)):
        eg.construct(tName, tStartsSets[i], blockSizesSets[i], exons)

    exonClusters, clusterReferences = eg.cluster(exons)
    geneModels = eg.buildGeneModels(exons, exonClusters, clusterReferences)
    starts, sizes = eg.printBed(geneModels)
    print >> sys.stderr, 'sizes = ', sizes
    print >> sys.stderr, 'starts = ', starts
    #assert sizes == ['15', '100', '50', '120', '50', '100']
    #assert starts == ['0', '55', '255', '455','655', '755']
    assert sizes == ['10', '50', '50', '120', '200']
    assert starts == ['0', '100', '250', '450', '650']

def testClusters3():
    exons = {}
    tStartsSets = [#[100, 300, 500, 700],
                    #[150, 300,500, 700,800],
                    [50, 150, 300, 700],
                    [50, 300, 500],
                    #[45, 100, 300, 500, 700, 800],
                    ]
    blockSizesSets = [#[100, 50, 120, 200],
                        #[50, 50, 150, 50, 100],
                        [10, 50, 320, 200],
                        [150, 50, 50],
                        #[15, 100, 50, 120, 50, 100],
                        ]
    tName = 'chr1'
    for i in range(len(tStartsSets)):
        eg.construct(tName, tStartsSets[i], blockSizesSets[i], exons)

    exonClusters, clusterReferences = eg.cluster(exons)
    geneModels = eg.buildGeneModels(exons, exonClusters, clusterReferences)
    starts, sizes = eg.printBed(geneModels)
    print >> sys.stderr, 'sizes = ', sizes
    print >> sys.stderr, 'starts = ', starts
    #assert sizes == ['15', '100', '50', '120', '50', '100']
    #assert starts == ['0', '55', '255', '455','655', '755']
    assert sizes == ['10', '50', '50', '120', '200']
    assert starts == ['0', '100', '250', '450', '650']

def testClusters4():
    exons = {}
    tStartsSets = [#[100, 300, 500, 700],
                    #[150, 300,500, 700,800],
                    #[50, 150, 300, 700],
                    [50, 300, 500],
                    [45, 100, 300, 500, 700, 800],
                    ]
    blockSizesSets = [#[100, 50, 120, 200],
                        #[50, 50, 150, 50, 100],
                        #[10, 50, 320, 200],
                        [150, 50, 50],
                        [15, 100, 50, 120, 50, 100],
                        ]
    tName = 'chr1'
    for i in range(len(tStartsSets)):
        eg.construct(tName, tStartsSets[i], blockSizesSets[i], exons)

    exonClusters, clusterReferences = eg.cluster(exons)
    geneModels = eg.buildGeneModels(exons, exonClusters, clusterReferences)
    starts, sizes = eg.printBed(geneModels)
    print >> sys.stderr, 'sizes = ', sizes
    print >> sys.stderr, 'starts = ', starts
    assert sizes == ['15', '100', '50', '120', '50', '100']
    assert starts == ['0', '55', '255', '455', '655', '755']

def testClusters5():
    exons = {}
    tStartsSets = [[100, 300, 500, 700],
                    #[150, 300,500, 700,800],
                    #[50, 150, 300, 700],
                    [50, 300, 500],
                    #[45, 100, 300, 500, 700, 800],
                    ]
    blockSizesSets = [[100, 50, 120, 200],
                        #[50, 50, 150, 50, 100],
                        #[10, 50, 320, 200],
                        [150, 50, 50],
                        #[15, 100, 50, 120, 50, 100],
                        ]
    tName = 'chr1'
    for i in range(len(tStartsSets)):
        eg.construct(tName, tStartsSets[i], blockSizesSets[i], exons)

    exonClusters, clusterReferences = eg.cluster(exons)
    geneModels = eg.buildGeneModels(exons, exonClusters, clusterReferences)
    starts, sizes = eg.printBed(geneModels)
    print >> sys.stderr, 'sizes = ', sizes
    print >> sys.stderr, 'starts = ', starts
    assert sizes == ['150', '50', '120', '200']
    assert starts == ['0', '250', '450', '650']

def testClusters6():
    exons = {}
    tStartsSets = [[100, 300, 500, 700],
                    [150, 300,500],
                    [100, 300, 700],
                    ]
    blockSizesSets = [[100, 50, 120, 200],
                        [50, 50, 150],
                        [100, 220, 200],
                        ]
    tName = 'chr1'
    for i in range(len(tStartsSets)):
        eg.construct(tName, tStartsSets[i], blockSizesSets[i], exons)

    exonClusters, clusterReferences = eg.cluster(exons)
    geneModels = eg.buildGeneModels(exons, exonClusters, clusterReferences)
    starts, sizes = eg.printBed(geneModels)
    print >> sys.stderr, 'sizes = ', sizes
    print >> sys.stderr, 'starts = ', starts
    assert sizes == ['100', '50', '120', '200']
    assert starts == ['0', '200', '400', '600']
