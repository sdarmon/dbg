import sys

Arg = sys.argv[:]
seq = []
vu = dict()
compt = 0
step = 0
with open(Arg[1], 'r') as f:
    for line in f:
        step = (step+1)%4
        if step == 0:
            if vu.get(line,-1) != -1:
                vu[line]+=1
            else:
                vu[line] = 1
            compt+=1
print("Il y a "+str(len(vu))+" under paths uniques sur les "+str(compt)+" bulles.")
with open(Arg[2], 'w') as o:
    for key,value in vu.items():
        o.write(key[:-1]+"\t"+str(value)+'\n')
