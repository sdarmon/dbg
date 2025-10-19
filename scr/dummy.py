import sys

Arg = sys.argv[:]

if len(Arg) != 4:
    print("Use : " + Arg[0] + "fileA fileB tableau")
    exit()

#Read fileA

id_A= set()
id_B= set()
with open(Arg[1], 'r') as f:
    for line in f:
        L = line[:-1]
        id_A.add(L)

#Read fileB
with open(Arg[2], 'r') as f:
    for line in f:
        L = line[:-1]
        id_B.add(L)

with open(Arg[3], 'r') as f:
    for line in f:
        L = line[:-1].split(',')
        seen_A = False
        seen_B = False
        for id in L:
            if id in id_A:
                seen_A = True
            if id in id_B:
                seen_B = True
        if seen_A and seen_B:
            print(line)