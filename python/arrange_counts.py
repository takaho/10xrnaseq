import os, sys, re, subprocess, argparse
import pandas as pd
import numpy as np

parser = argparse.ArgumentParser()
parser.add_argument('-i', nargs='+', metavar='directory', help='SAM file containers')
parser.add_argument('-p', default='/home/takaho/share/build/subread-1.6.0-source/bin/featureCounts', help='featureCount path')
parser.add_argument('-g', default='/home/takaho/share/ngs/mm10/genes.gtf', metavar='GTF file', help='annotation')
parser.add_argument('-o', metavar='direcotory', help='output directory')
parser.add_argument('--batch_size', type=int, default=100, help='featureCount input size')
parser.add_argument('--num_batches', type=int, default=0, help='How many batches will be used')
parser.add_argument('--log', default='log.txt')

args = parser.parse_args()
prog = args.p#'/home/takaho/share/build/subread-1.6.0-source/bin/featureCounts'
gtf = args.g#'/home/takaho/share/ngs/mm10/genes.gtf'
dstdir = args.o#'test_results'
integration = './cxx/normalize_fc'
#limit = 1000
batch_size = args.batch_size
num_batches = args.num_batches
if os.path.exists(dstdir) is False: os.makedirs(dstdir)

log = open(args.log, 'w')
for dn in args.i:
    filenames = []
    fn_out = []
    bn = os.path.basename(dn).strip('_')
    fn_out = os.path.join(dstdir, bn + '.fc')
    log.write(dn + '\n')
    if not os.path.exists(fn_out) or os.path.getsize(fn_out) <= 10000:
        log.write('exec featureCount\n')
        for fn in os.listdir(dn):
            if fn.endswith('.sam'):
                filenames.append(os.path.join(dn, fn))
        start = 0
        slots = []
        num_slots = len(filenames) // batch_size
        while start < len(filenames):
            stop = min(len(filenames), start + batch_size)
            slots.append(filenames[start:stop])
            start += batch_size
            if num_batches > 0 and len(slots) >= num_batches: break
        output_filenames = []
        for i, fns in enumerate(slots):
            output_filename = os.path.join(dstdir, bn + '.{}.fc'.format(i))
            cmd = [prog, '-T', '4', '-o', output_filename, '-a', gtf ] + fns#filenames
            log.write('{}\t{}\n'.format(output_filename, len(fns)))
            subprocess.Popen(cmd).wait()
            if os.path.exists(output_filename) and os.path.getsize(output_filename) > 10000:
                output_filenames.append(output_filename)
            else:
                log.write('ERROR:{}\n'.format(' '.join(cmd)))
        log.write('merge {} batch files\n'.format(len(filenames)))
        log.flush()

        if 1:
            cmd = [integration, '-o', os.path.join(dstdir, bn), '-i', ] + output_filenames
            log.write(' '.join(cmd) + '\n')
            sys.stderr.write(' '.join(cmd) + '\n')
            log.flush()
            subprocess.Popen(cmd).wait()
            continue
                
        header_items = []
        table = []
        for fn in output_filenames:
            if os.path.exists(fn) and os.path.getsize(fn) > 100000:
                with open(fn) as fi:
                    while 1:
                        line = fi.readline()
                        if line.startswith('#'):
                            continue
                        break
                    helems = [x_.split('/')[-1].split('.')[0] for x_ in line.strip().split('\t')]
                    if len(header_items) == 0:
                        header_items = helems
                        cols_start = 0
                    else:
                        cols_start = 6
                        header_items += helems[cols_start:]
                    row = 0
                    for line in fi:
                        items = line.strip().split('\t')
                        if cols_start == 0:
                            items[1] = items[1].split(';', 1)[0]
                            start = [int(x_) for x_ in items[2].split(';')]
                            end = [int(x_) for x_ in items[3].split(';')]
                            items[2] = '{}'.format(min(start))
                            items[3] = '{}'.format(max(end))
                            items[4] = items[4][0] # ori
                            table.append(items)
                        else:
                            table[row] += items[cols_start:]
                        row += 1
#                for f_ in (fn, fn + '.summary'):
#                    os.unlink(f_)
        with open(fn_out, 'w') as fo:
            fo.write('#featureCount of SAM files in {}, with {}\n'.format(dn, gtf))
            fo.write('\t'.join(header_items) + '\n')
            for row in table:
                fo.write('\t'.join(row) + '\n')
                    

    label = os.path.basename(dn).strip('_')
    m = re.search('lib\\d+', label)
    if m:
        label = m.group(0)
    if os.path.exists(fn_out) and os.path.getsize(fn_out) > 10000:
        table = pd.read_csv(fn_out, sep='\t', skiprows=1, index_col=0)
        cnt = table[table.columns[6:]].copy()
        cols = ['{}_{}'.format(label, i + 1) for i in range(len(cnt.columns))]
        cnt.columns = cols
        fn_cnt = os.path.join(dstdir, os.path.basename(dn).strip('_') + '.cnt')
        fn_tpm = os.path.join(dstdir, os.path.basename(dn).strip('_') + '.tpm')
        cnt.to_csv(fn_cnt, sep='\t')
        for col in cols:
            cnt[col] *= 1e6 / np.sum(cnt[col])
        cnt.to_csv(fn_tpm, sep='\t')
    log.flush()
        
log.close()
