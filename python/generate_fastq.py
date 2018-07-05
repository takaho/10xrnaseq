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

              
import numpy as np
for name in dataset.keys():
    dstdir = 'out_{}'.format(name)
    if os.path.exists(dstdir) is False: os.makedirs(dstdir)
    for fn in os.listdir(dstdir):
        os.unlink(os.path.join(dstdir, fn))
    
    filestreams = {}
    procs = []
    istrs = []
    for section in ('R1', 'R2'):
#    for section in ('I1', 'R1', 'R2'):
        cmd = 'zcat', dataset[name][section]
        proc = subprocess.Popen(cmd, stdout=subprocess.PIPE)
        procs.append(proc)
        istrs.append(proc.stdout)
#    print(name)
#    istrs = [gzip.open(dataset[section]) for section in ('I1', 'R1', 'R2')]
    prev = ''
    num_reads = 0
    freq = {}
    def flush_and_drop(filestreams, tags=None):
        if tags is None:
            for fh, ln in filestreams.values():
                for l in ln:
                    fh.write(l)
                fh.close()
            for t in filestream.keys():
                filestreams.drop(t)
                pass
        else:
            for t in tags:
                fh, ln = filestreams[t]
                for l in ln:
                    fh.write(l)
                fh.close()
                filestreams.pop(t)
#    max_handlers = cat /proc/sys/fs/file-max
    #max_handlers = 3252599
    max_handlers = 1000#4096#16384
    #128
    while 1:
        data = []
        for i in range(2):
            lines = []
            for j in range(4):
                lines.append(istrs[i].readline().decode('utf-8'))
            data.append(lines)
        if data[0][0] == '': break
        tag = data[0][0].split(':')[-1].strip() + '_' + data[0][1].strip()
        if tag not in filestreams:
            N = len(filestreams)
            if N >= max_handlers:
                sys.stderr.write('{} {} => '.format(num_reads, N))
                ts = list(filestreams.keys())
                closing = sorted(filestreams.keys(), key=lambda k:len(filestreams[k][1]), reverse=True)[0:N-max_handlers // 2]
                #closing = set([ts[np.random.randint(0, N)] for i in range(N - max_handlers // 2)])
                flush_and_drop(filestreams, closing)
                sys.stderr.write('{}\n'.format(len(filestreams)))
            fn_out = os.path.join(dstdir, tag + '.fastq')
            filestreams[tag] = [open(fn_out, 'a'), data[1]]
        else:
            filestreams[tag][1] += data[1]
            pass
        num_reads += 1
    flush_and_drop(filestreams)
    
#        freq[tag] = freq.get(tag, 0) + 1
        
#         if prev != tag:
#             freq[tag] = freq.get(tag, 0) + 1
#             tags = [data[i][0].split(':')[-1].strip() for i in range(3)]
#             if tag != tags[0] or tag != tags[1] or tag != tags[2]:
#                 print('{}, {}/{}/{}\t{}'.format(tag, tag == tags[0], tag == tags[1], tag == tags[2], repr(tags)))
# #            print('{} {}, {}, {}, {}'.format(num_reads, data[0][0].strip(), data[1][0].strip(), data[2][0].strip(), tag))
#             prev = tag
#         else:
#             freq[tag] += 1

    #     if num_reads % 1000000 == 0:
    #         for t in sorted(freq.keys(), key=lambda t:freq[t], reverse=True):
    #             n = freq[t]
    #             if n > 100:
    #                 print('{}\t{}\t{:.3f}%'.format(t, n, n * 100.0 / num_reads))
            
    # for istr in istrs:
    #     istr.close()
    # with open('{}.r1'.format(name), 'w') as fo:
    #     for t in sorted(freq.keys(), key=lambda t:freq[t], reverse=True):
    #         n = freq[t]
    #         fo.write('{}\t{}\n'.format(t, n))

        
            

