#!/usr/bin/env python3
import sys

if len(sys.argv) > 1:
    infile = open(sys.argv[1],'r')
else:
    print('Specify FASTA input file.')

lines = infile.readlines()

for i in range(0,len(lines)):
    if lines[i].startswith('>'):
        print(lines[i].strip())
    elif i < len(lines)-1 and lines[i+1].startswith('>'):
        print(lines[i].strip())
    else:
        print(lines[i].strip(),end="", flush=True)
