import os, sys, re, argparse, gzip

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


parser = argparse.ArgumentParser()
parser.add_argument('-b', nargs='+')
parser.add_argument('--R1', nargs='+')
parser.add_argument('--R2', nargs='+')
parser.add_argument('-o', default='out')
parser.add_argument('--minimum-reads', type=int, default=1000)
parser.add_argument('--verbose', action='store_true')
args = parser.parse_args()

size = args.minimum_reads
dstdir = args.o
verbose = args.verbose
if verbose:
    lap_period = 1000000
else:
    lap_period = 0
#    lap_period = 1
if os.path.exists(dstdir) is False:
    os.makedirs(dstdir)
    
# get whitelist
counts = {}
for fn in args.b:
    with open(fn) as fi:
        for line in fi:
            items = line.strip().split('\t')
            if line.startswith('#') is False and len(items) >= 2 and items[1].isdigit():
                counts[items[0]] = counts.get(items[0], 0) + int(items[1])

applicable = set([b for b in sorted(counts.keys(), key=lambda b:counts[b], reverse=True) if counts[b] >= size])
second_code = {}
for c in applicable:
    second_code[get_code(c)[0]] = c
    
print(len(applicable))
#exit()
#outputfiles = []
N = len(args.R1)
if N != len(args.R2):
    raise Exception('not balance')

indexes = """
ACAGAGGT TATAGTTG CGGTCCCA GTCCTAAC
GCATCTCC TGTAAGGT CTGCGATG AACGTCAA
GGTTTACT CTAAACGG TCGGCGTC AACCGTAA
CAGTACTG AGTAGTCT GCAGTAGA TTCCCGAC
"""
well_indexes = re.split('[^ACGT]+', indexes.strip())
wellnorm = {}
for i, c in enumerate(well_indexes):
    bincode = get_code(c)[0]
    wellnorm[bincode] = c

wells = {}
single_cells = {}
cache_size = 1000
handlers = []
seq_cache = {}
from gmpy2 import popcount
for i in range(N):
    z1 = gzip.open(args.R1[i])
    z2 = gzip.open(args.R2[i])
    num_reads = 0
    while 1:
        l1 = [z1.readline() for i in range(4)]
        l2 = [z2.readline() for i in range(4)]
        num_reads += 1
        if lap_period > 0 and num_reads % lap_period == 0:
            print('{} // {}'.format(num_reads // lap_period, len(single_cells)))
            for key in sorted(single_cells.keys(), key=lambda s:single_cells[s], reverse=True)[0:20]:
                print('{}\t{}'.format(key, single_cells[key]))
        if l1[-1] == '': break
        h = l1[0].decode('utf-8')

        well_code = h[h.rfind(':') + 1:].strip()
        if 1:
            bincode, binmask = get_code(well_code)
            num_N = popcount(binmask) >> 1
            well_code = None
            for bc, c in wellnorm.items():
                op_ = (bc ^ bincode) & 0xffff
                num_mm = popcount(((op_ & 0x5555) | ((op_ >> 1) & 0x5555))) + num_N
                if num_mm <= 1:
                    well_code = c
                    break
            if well_code is None:
                continue
            pass
        elif well_code not in well_indexes:
            continue
        barcode = l1[1].decode('utf-8').strip()
        cell_code = barcode[0:16]
        umi = barcode[16:]
        if 0:
            bincode, binmask = get_code(cell_code)
            detected = None
            num_N = popcount(binmask) >> 1
            for bc, c in second_code.items():
                op_ = (bc ^ bincode) & 0xffffffff
                num_mm = popcount(((op_ & 0x55555555) | ((op_ >> 1) & 0x55555555))) + num_N
                if num_mm <= 1:
                    detected = c
                    break
            if detected is None:
                code_ = well_code + '_?'
                single_cells[code_] = single_cells.get(code_, 0) + 1
                #print(well_code)
                continue
            cell_code = detected
        if cell_code in applicable:
            wells[cell_code] = wells.get(cell_code, 0) + 1
            code_ = well_code + '_' + cell_code
            if code_ not in single_cells:
                single_cells[code_] = 0
                filename = os.path.join(dstdir, code_ + '.fastq')
                open(filename, 'w').close()
            else:
                single_cells[code_] += 1#= single_cells.get(code_, 0) + 1
            if len(seq_cache) >= cache_size:
                for c_, seqs in seq_cache.items():
                    fn = os.path.join(dstdir, c_ + '.fastq')
                    with open(fn, 'a') as fo:
                        start = single_cells[c_] - len(seqs) + 2
                        for i, seq in enumerate(seqs):
                            fo.write('@{}_{}:'.format(c_, start + i))
                            fo.write(seq)#+ seq)
                            pass
                        pass
                seq_cache = {}
            if code_ not in seq_cache:
                seq_cache[code_] = []
            seq_cache[code_].append(umi + '\n' + l2[1].decode('utf-8')[1:] + '+\n' + l2[3].decode('utf-8'))
        else:
             code_ = well_code + '_?'
             single_cells[code_] = single_cells.get(code_, 0) + 1
        pass
    z1.close()
    z2.close()
    for c_, seqs in seq_cache.items():
        fn = os.path.join(dstdir, c_ + '.fastq')
        with open(fn, 'a') as fo:
            start = single_cells[c_] - len(seqs)
            for i, seq in enumerate(seqs):
                fo.write('@{}_{}:'.format(c_, start + i))
                fo.write(seq)#+ seq)
                pass
            pass
