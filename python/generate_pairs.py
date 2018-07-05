import os, sys, re, subprocess

srcdirs = '1804KHF-0105_10x_RawData_Outs', '1804KHF-0146_10x_RawData_Outs'

minimum = 1000
mismatches = 1
prog = './cxx/a.out'
for srcdir in srcdirs:
    for path, dns, fns in os.walk(srcdir):
        for fn in fns:
            m = re.match('(lib\\d+)_.*_R1_.*\\.fastq\\.gz$', fn)
            if m:
                #print(m)
                r1 = os.path.join(path, fn)
                #print(r1)
                r2 = os.path.join(path, fn.replace('_R1_', '_R2_'))
                if os.path.exists(r2):
                    #print(r2)
                    mh = re.search('\\W(HM\\w+)', path)
                    if mh:
                        run = mh.group(1)
                        lib = m.group(1)
                        #print(run, lib)
                        fn_out = 'm{}_k{}_{}_{}.fastq.gz'.format(minimum, mismatches, run, lib, run)
                        cmd = prog, '--R1', r1, '--R2', r2, '-m', '1000', '-k', '{}'.format(mismatches), '-o', fn_out, '--verbose', '-L', '10000'
                        #print(' '.join(cmd))
                        subprocess.Popen(cmd).wait()
