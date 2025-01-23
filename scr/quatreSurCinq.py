import sys

Arg = sys.argv[:]
step = 0
with open(Arg[1], 'r') as f:
    for line in f:
        step = (step+1)%5
        if step != 0:
            print(line[:-1])
