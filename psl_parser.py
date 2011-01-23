#! /usr/local/bin/python
'''
    This parser is used to parse PSL file (i.e. from BLAT output).
    read() return an object containing all information of one line of PSL format.
'''

import sys

class PSL(object):
    def __init__(self, **kwargs):
        self.attrib = kwargs

def read(fobj, comment):
    '''
        fobj = file object.
        comment = read() will ignore the line starting with comment character.
    '''
    for line in fobj:
        if line.startswith(comment): continue
        attrib = {}
        rows = line.split()
        attrib['matches'] = rows[0]
        attrib['misMatches'] = rows[1]
        attrib['repMatches'] = rows[2]
        attrib['nCount'] = rows[3]
        attrib['qNumInsert'] = rows[4]
        attrib['qBaseInsert'] = rows[5]
        attrib['tNumInsert'] = rows[6]
        attrib['tBaseInsert'] = rows[7]
        attrib['strand'] = rows[8]
        attrib['qName'] = rows[9]
        attrib['qSize'] = rows[10]
        attrib['qStart'] = rows[11]
        attrib['qEnd'] = rows[12]
        attrib['tName'] = rows[13]
        attrib['tSize'] = rows[14]
        attrib['tStart'] = rows[15]
        attrib['tEnd'] = rows[16]
        attrib['blockCount'] = rows[17]
        attrib['blockSizes'] = rows[18].split(',')[:-1]
        attrib['qStarts'] = rows[19].split(',')[:-1]
        attrib['tStarts'] = rows[20].split(',')[:-1]
        
        pobj = PSL(**attrib) # pobj = PSL object
        yield pobj

if __name__ == '__main__':
    pass
