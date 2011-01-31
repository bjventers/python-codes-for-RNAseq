import sys

fp1 = sys.argv[1]
fp2 = sys.argv[2]
snp1 = sys.argv[3]
snp2 = sys.argv[4]
s1 = {}
s2 = {}
for line in open(snp1):
    chr, pos, ref, geno, cov = line.strip().split()
    pos = int(pos)
    s1[pos] = geno

for line in open(snp2):
    chr, pos, ref, geno, cov = line.strip().split()
    pos = int(pos)
    s2[pos] = geno
    
juncs = {}
for line in open(fp1):
    rows = [int(l) for l in line.strip().split(',')[:-1]]
    start = rows[0]
    juncs[start] = rows[1:]

for line in open(fp2):
    rows = [int(l) for l in line.strip().split(',')[:-1]]
    start = rows[0]
    try:
        j1 = set(juncs[start])
        j2 = set(rows[1:])
        diff1 = j1.difference(j2) 
        diff2 = j2.difference(j1)
        if diff1 or diff2:
            print start,
            if diff1: 
                print '\t',
                for s in diff1:
                    if s in juncs[start]: subset = 1
                    else: subset = 2
                    print '%d:%d' % (subset,s),
                    for j in range(s-3, s+3):
                        if j in s1 and j in s2: continue
                        elif j in s1:
                            print '1:%d:%s;' % (j, s1[j]),
                        elif j in s2:
                            print '2:%d:%s;' % (j, s2[j]),
                    print ',',
                print '\t',
            else: print '\tNone',
            if diff2: 
                print '\t',
                for s in diff2:
                    if s in juncs[start]: subset = 1
                    else: subset = 2
                    print '%d:%d' % (subset,s),
                    for j in range(s-3, s+3):
                        if j in s1 and j in s2: continue
                        elif j in s1:
                            print '1:%d:%s;' % (j, s1[j]),
                        elif j in s2:
                            print '2:%d:%s;' % (j, s2[j]),
                    print ',',
            else: print '\tNone',
            print
    except KeyError:
        pass
