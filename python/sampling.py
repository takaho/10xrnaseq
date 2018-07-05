import os, sys, re
import subprocess, math
import argparse
import numpy as np

parser = argparse.ArgumentParser()
parser.add_argument('-i', metavar='sam file')
parser.add_argument('-n', type=int, default=1000000)
parser.add_argument('--m', type=int, default=1)
parser.add_argument('--batch', type=int, default=100)
parser.add_argument('--batch-param', type=int, default=10)
args = parser.parse_args()


fn = args.i
size = os.path.getsize(args.i)
batch = args.batch
bparam = args.batch_param

fh = open(fn)
while 1:
    pos = np.random(0, size - batch * 256)
    fh.seek(pos)
    for i in range(batch):
        l = fh.readline()
        if l.startswith('@'):
        if np.random.randint(0, bparam) == 0:
            
    
    
    
num
fh = open(sys.argv[1])

with open(sys.argv[1]) as fi:
