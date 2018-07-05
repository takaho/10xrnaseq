import os, sys, re, subprocess, pathlib, gzip

"""1804KHF-0105_10x_RawData_Outs/fastq_path/HMMWYBBXX/lib7/lib7_S1_L002_I1_001.fastq.gz"""

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

        
            

