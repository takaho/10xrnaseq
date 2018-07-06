import os, sys, re, argparse, subprocess, pathlib

"""
Procedure

ount_and_rename --R1 [R1.fastq.gz] --R2 [R2.fastq.gz] -o [renamed_fastq.gz] --mismatches [num] 
=> fastq.gz
STAR 
=> BAM
python split_sam.py -i [BAM] -o [output_directory] -n [minimum reads]
=> many sam files
featureCount
=> .fc files
normalize_fc -i [fc2 files] -o [label]
=> label.cnt and label.tpm files


count_and

"""

parser = argparse.ArgumentParser()
parser.add_argument('-i', metavar='bam/sam file')
parser.add_argument('-o', metavar='directory')
parser.add_argument('-n', type=int, default=5000, metavar='number', help='minimum read counts for a cell')
parser.add_argument('--forced', action='store_true')
parser.add_argument('--verbose', action='store_true')
args = parser.parse_args()

filename_sam = args.i
if filename_sam.lower().endswith('.bam'):
    cmd = 'samtools', 'view', '-h', filename_sam
    proc = subprocess.Popen(cmd, stdout=subprocess.PIPE)
    istr = proc.stdout
else:
    proc = None
    istr = open(filename_sam)
    pass

umi_cache = {}
current_pos = -1
header = ''
dstdir = pathlib.Path(args.o)
seqs = {}
cache_size = 10000000
forced = args.forced
verbose = args.verbose
sample_read_count = {}
num_reads = 0            
num_cached_reads = 0
lap = 100000
minimum_counts = args.n#reads = args.n
num_duplicated = 0

def dump_sequences(seqs, dstdir, header, threshold=500, check_file=False, verbose=False):
    num_saved = 0
    num_kept = 0
    for key, val in seqs.items():
        mode = 0
        num = len(val)
        if num == 0:
            mode = 0
        elif check_file:
            dstfile = dstdir.joinpath(key + '.sam')
            if dstfile.exists():
                mode = 2
            elif num >= threshold:
                mode = 1
        else:
            if num >= threshold:
                dstfile = dstdir.joinpath(key + '.sam')
                if dstfile.exists():
                    mode = 2
                else:
                    mode = 1
        if mode == 1:
            fh = dstfile.open('w')
            fh.write(header)
        elif mode == 2:
            fh = dstfile.open('a')
        else:
            num_kept += num
            continue
        num_saved += num
        for s in val: fh.write(s)
        fh.close()
        seqs[key] = []
    if verbose:
        sys.stderr.write(' saved:{} / kept:{}       \r'.format(num_saved, num_kept))
        
#minimum_reads = 500                #for l in header

if dstdir.exists():
    if not forced:
        raise Exception('directory {} exists'.format(dstdir.as_posix()))
    if verbose:
        sys.stderr.write(' '.join(['rm', '-rf', dstdir.as_posix()]) + '\n')
    subprocess.Popen(['rm', '-rf', dstdir.as_posix()]).wait()
    #if not dstdir.exists():
    pass
dstdir.mkdir(exist_ok=True, parents=True)

while 1:
    line = istr.readline()
    if line == '': break
    line = line.decode('utf-8')
    if line.startswith('@'):
        header += line
        continue
    items = line.split('\t')
    if len(items) < 4:
        sys.stderr.write('terminate by ' + line)
        break
    if items[2] == '*' or items[2].find('_') >= 0:
        if proc is not None:
            proc.terminate()
#        istr.close()
        continue
    num_reads += 1
    #num_cached_reads += 1
    if verbose:
        if lap > 0 and num_reads % lap == 0:
            sys.stderr.write(' {} {} {}       \r'.format(num_reads // lap, items[2], items[3]))
    pos = int(items[3])
    try:
        tag, num, umi = items[0].split(':')
    except:
        continue
    duplicated = False
    if pos == current_pos:
        if tag in umicache:
            if umi in umicache[tag]: # duplicated
                #sys.stderr.write('{} duplicated       \n'.format(items[0]))
                num_duplicated += 1
                duplicated = True
            else:
                umicache[tag].append(umi)
        else:
            umicache[tag] = [umi,]
    else:
        umicache = {}
        current_pos = pos
    if not duplicated:
        if num_cached_reads >= cache_size:
            #if verbose:
            #    sys.stderr.write('dump {} sequences  \r'.format(num_cached_reads))
            dump_sequences(seqs, dstdir, header, minimum_counts // 100)
            num_cached_reads = sum([len(x_) for x_ in seqs.values()])
        if tag not in seqs:
            seqs[tag] = [line,]
            sample_read_count[tag] = 1
        else:
            seqs[tag].append(line)
            sample_read_count[tag] += 1
        num_cached_reads += 1
if verbose: sys.stderr.write('complete processes               \n')
accepted = []
seqs_to_save = {}
num_cells = 0
for key, val in sample_read_count.items():
    if val >= minimum_counts:
        seqs_to_save[key] = seqs[key]
        num_cells += 1
dump_sequences(seqs_to_save, dstdir, header, 0)
istr.close()
if proc is not None:
    proc.wait()
    
num_accepted = num_rejected = 0
num_removed = 0
for sample, n in sample_read_count.items():
    if n < minimum_counts:
        num_rejected += n
        dstfile = dstdir.joinpath(sample + '.sam')
        num_removed += 1
        if dstfile.exists():
            dstfile.unlink()
    else:
        num_accepted += n
        
if verbose:
    sys.stderr.write('filename:{}\noutput:{}\nremoved_files:{}\nsaved cells:{}\nreads accepted:{}\nreads rejected:{}\nreads duplicated:{}\nreads total:{}\n'
                     .format(filename_sam, dstdir.as_posix(), num_removed, num_cells, num_accepted, num_rejected, num_duplicated, num_reads))
    

