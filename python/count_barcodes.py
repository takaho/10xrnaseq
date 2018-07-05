import os, sys, re, subprocess, pathlib, gzip

"""1804KHF-0105_10x_RawData_Outs/fastq_path/HMMWYBBXX/lib7/lib7_S1_L002_I1_001.fastq.gz"""


first_indices = {
    'AACGTCAA':0,#	85457155
    'GCATCTCC':1,#	84432940
    'CTGCGATG':2,#	70279153
    'TGTAAGGT':3,#	54412965

    'GTCCTAAC':4,#	83802878
    'CGGTCCCA':5,#	78867120
    'ACAGAGGT':6,#	73400348
    'TATAGTTG':7}#	67410698

r1 = {}
for path, dns, fns in os.walk('1804KHF-0105_10x_RawData_Outs'):
    for fn in fns:
        m = re.match('(lib\\d+_S\\d+).*_(I1|R1|R2)_\\d{3}\\.fastq.gz$', fn)
        if m is not None and m.group(2) == 'R1':
            r1[m.group(1)] = os.path.join(path, fn)
            #:
            # print(m.groups(), fn)
            # name = m.group(1)
            # if name not in dataset:
            #     dataset[name] = {}#[None, None, None]
            # dataset[name][m.group(2)] = os.path.join(path, fn)
num_reads = 0
for name, fn in r1.items():
    freq = {}
    with gzip.open(fn) as fi:
        while 1:
            n = fi.readline()
            if n == '': break
            s = fi.readline().decode('utf-8')
            fi.readline()
            fi.readline()
            code = n[-9:-1].decode('utf-8')
            idx1 = None
            for idx, num in first_indices.items():
                mm = 0
                if idx == code:
                    idx1 = idx
                    break
    #            print(idx, s)
                for i, c in enumerate(code):
    #                print(i, c, s, idx)
                    if c != idx[i]:
                        mm += 1
                        if mm >= 2:
                            break
                if mm < 2:
                    idx1 = idx
                    break
            if idx1 is not None:
                idx2 = s[0:-1]
                barcode = idx1 + '_' + idx2
                freq[barcode] = freq.get(barcode, 0) + 1
            num_reads += 1
            if num_reads % 100000 == 0:
                for bc in sorted(freq.keys(), key=lambda b:freq[b], reverse=True)[0:20]:
                    print('{}\t{}\t{:.4f}'.format(bc, freq[bc], freq[bc] * 100.0 / num_reads))
                
with open(name + '.freqall', 'w') as fo:
    for bc in sorted(freq.keys(), key=lambda b:freq[b], reverse=True)[0:20]:
        fo.write('{}\t{}\t{:.4f}\n'.format(bc, freq[bc], freq[bc] * 100.0 / num_reads))
            
        
exit()

dataset = {}
for path, dns, fns in os.walk('.'):
    for fn in fns:
        m = re.match('(lib\\d+_S\\d+).*_(I1|R1|R2)_\\d{3}\\.fastq.gz$', fn)
        if m:
            print(m.groups(), fn)
            name = m.group(1)
            if name not in dataset:
                dataset[name] = {}#[None, None, None]
            dataset[name][m.group(2)] = os.path.join(path, fn)
print(dataset)
for name in dataset.keys():
    procs = []
    istrs = []
    for section in ('I1', 'R1', 'R2'):
        cmd = 'zcat', dataset[name][section]
        proc = subprocess.Popen(cmd, stdout=subprocess.PIPE)
        procs.append(proc)
        istrs.append(proc.stdout)
#    print(name)
#    istrs = [gzip.open(dataset[section]) for section in ('I1', 'R1', 'R2')]
    prev = ''
    num_reads = 0
    freq = {}
    while 1:
        data = []
        for i in range(3):
            lines = []
            for j in range(4):
                lines.append(istrs[i].readline().decode('utf-8'))
            data.append(lines)
        tag = data[0][1].strip()
        if tag == '':
            break
        num_reads += 1
        if prev != tag:
            freq[tag] = freq.get(tag, 0) + 1
            tags = [data[i][0].split(':')[-1].strip() for i in range(3)]
            if tag != tags[0] or tag != tags[1] or tag != tags[2]:
                print('{}, {}/{}/{}\t{}'.format(tag, tag == tags[0], tag == tags[1], tag == tags[2], repr(tags)))
#            print('{} {}, {}, {}, {}'.format(num_reads, data[0][0].strip(), data[1][0].strip(), data[2][0].strip(), tag))
            prev = tag
        else:
            freq[tag] += 1
        if num_reads % 1000000 == 0:
            for t in sorted(freq.keys(), key=lambda t:freq[t], reverse=True):
                n = freq[t]
                if n > 100:
                    print('{}\t{}\t{:.3f}%'.format(t, n, n * 100.0 / num_reads))
            
    for istr in istrs:
        istr.close()
    with open('{}.freq'.format(name), 'w') as fo:
        for t in sorted(freq.keys(), key=lambda t:freq[t], reverse=True):
            n = freq[t]
            fo.write('{}\t{}\n'.format(t, n))

        
            

