#!/usr/bin/env python3

import sys

if len(sys.argv)<2:
    print ("usage:", sys.argv[0], 'input > output')
    exit(0)
fname = sys.argv[1]
f = open(fname, 'r')
contents = f.read()
f.close()

contents = contents.split('\n')
contents = list(filter(lambda x: len(x) and x[0]!='!', contents))

qs = []
atom_start = []
for i in range(len(contents)):
    line = contents[i].split()
    if len(line)==2 and line[0]=='a' and line[1].isdigit():
        atom_start.append(i)
        qs.append(int(line[1]))

if atom_start[0]!=0:
    exit(1)
atom_start.append(len(contents))

basis = {}
for i in range(len(atom_start)-1):
    basis[qs[i]] = contents[ atom_start[i]+1 : atom_start[i+1] ]


print('$basis')
print('type=gc')
for q in basis.keys():
    nl = 0
    j = 0
    while j<len(basis[q]):
        line = basis[q][j].split()
        if len(line)!=3 and line[0]!='H' and not line[1].isdigit() and not line[2].isdigit():
            exit(1)
        nl += 1
        np = int(line[1])
        nc = int(line[2])
        basis[q][j] = '%2d %2d' % (np, nc)
        j+=1

        for i in range(np):
            line = basis[q][j].split()
            if len(line)<nc+1:
                exit(1)
            line = " ".join(['%+.12e'%float(x) for x in line])
            basis[q][j] = line
            j+=1

    print("      %3d %d" % (q, nl-1))
    print( "\n".join(basis[q]) )
print('$end')

