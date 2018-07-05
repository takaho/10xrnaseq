import os, sys, re
import numpy
import time
import gmpy2

def p1(n):
    return numpy.bincount(n)

def p2(n):
    return gmpy2.popcount(n)

def p3(n):
    return bin(n).count('1')

num_loops = 100
for i in range(num_loops):
    r = numpy.random.randint(0, 10000)
    for j in range(1000):
        c = p1(r)
    print(r, c)
    
    
