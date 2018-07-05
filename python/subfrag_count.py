import os, sys, re
import argparse
import gmpy2

def get_code(seq):
    num = 0
    mask = 0
    for i in range(len(seq)):
        s = seq[len(seq) - 1 - i]
#    for s in seq:
        n = m = 0
        if s == 'C':
            n = 1
        elif s == 'G':
            n = 2
        elif s == 'T':
            n = 3
        elif s == 'N':
            m = 3
        num = (num << 2) | n
        mask = (mask << 2) | m
    return num, mask
def get_mm(seq1, seq2):
    n = 0
    for i in range(len(seq1)):
        if seq1[i] != seq2[i]:
            n += 1
    return n
#def popcount()
#GGTTTACT
indexes = """
ACAGAGGT TATAGTTG CGGTCCCA GTCCTAAC
GCATCTCC TGTAAGGT CTGCGATG AACGTCAA
GGTTTACT CTAAACGG TCGGCGTC AACCGTAA
CAGTACTG AGTAGTCT GCAGTAGA TTCCCGAC
"""
wellindexes = ['AACGTCAA',
               'GCATCTCC',
               'CTGCGATG',
               'TGTAAGGT',

               'GTCCTAAC',
               'CGGTCCCA',
               'ACAGAGGT',
               'TATAGTTG']


indexes = re.split('[^ACGT]+', indexes.strip())
sequences = [get_code(seq)[0] for seq in indexes]
wellcodes = [get_code(seq)[0] for seq in wellindexes]

parser = argparse.ArgumentParser()
parser.add_argument('-i', metavar='filename')
parser.add_argument('-o', metavar='filename')
args = parser.parse_args()
fn = args.i
mismatches = 3

freq = [0] * (len(sequences) + 2)
num_reads = 0
lap = 10000#0000
mask = 0xffff
bmask = 0xaaaa

def popcount(n):
    return gmpy2.popcount(n)

def popcount1(n):
    return gmpy2.popcount(n)

def popcount2(n):
    return bin(n).count('1')
#    c = (n & 0x5555) + ((n >> 1) & 0x5555)
#    c = (c & 0x3333) + ((c >> 2) & 0x3333)
#    c = (x & 0x0f0f) + ((
wfreq = [0] * (len(wellcodes) + 1)
with open(fn) as fi:
    while 1:
        lines = [fi.readline() for i in range(4)]
        if lines[0] == '': break
        welltag = lines[0].split(':')[-1].strip()
        wcode, mask = get_code(welltag)
        num_N = popcount(mask) >> 1
        bmask = ~mask
        wellid = -1
        for i, w in enumerate(wellcodes):
            num_mm = popcount((wcode ^ w) & bmask) + num_N
            if num_mm < 2:
                wellid = i
                break
        wfreq[wellid] += 1
        #print(wellid, welltag, wellindexes[wellid])
        
        seq = lines[1].strip()
        N = len(seq)
        code, npat = get_code(seq)
        found = [[0] * len(indexes) for x_ in range(mismatches + 1)]
        for i in range(N - 8):
            m = (npat >> (i << 1)) & mask # N
            n = ~m
            s = (code >> (i << 1)) & mask # non N
            #pm = popcount(n)
            #num_N = gmpy2.popcount(m) >> 1
            num_N = popcount1(m) >> 1
#            frag = seq[i:i+8]
            #print('{} {:x} {:x} {:x} {:x}'.format(seq, code, m, n, s))
            for j, barcode in enumerate(sequences):
                a = (barcode ^ s) & n
                num_mm = popcount1((a | (a << 1)) & bmask) + num_N
                #num_mm = gmpy2.popcount((a | (a << 1)) & bmask) + num_N
                if num_mm <= mismatches:
                    found[num_mm][j] = 1
        best = -1 # not found
        for i, f in enumerate(found):
#            print(i, f)
            for j, v in enumerate(f):
                if v > 0:
                    if best != -1:
                        best = -2
                        break
                    else:
                        best = j
            if best != -1:
                break
#        print(best)
        freq[best] += 1
                        
        num_reads += 1
        if lap > 0 and num_reads % lap == 0:
            print(' {} : {}'.format(num_reads // lap, ' '.join(['{}'.format(x) for x in freq])))
            print(wfreq)
            if num_reads > 10 * lap: exit()
        
        
