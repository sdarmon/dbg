# Petit script permettant de supprimer les lignes ayant les mêmes target
# pour intersectionTE et intersectionKiss

import numpy as np 
import sys


Arg = sys.argv[:]

if len(Arg) not in [3,5]:
    print("Use : "+Arg[0]+ " input.txt output.txt -s pos")
elif len(Arg) == 5 and Arg[3] == "-s": #Cas Kissplice :
    with open(Arg[1],'r') as f:
        with open(Arg[2],'w') as o:
            Vu = []
            for line in f:
                if line[0] == "#":
                    continue
                target = int((line.split("\t")[int(Arg[4])]).split("_")[1]) #On récupère le numéro de la séquence
                if target >= len(Vu):
                    for i in range(len(Vu),target+2): #On cherche à savoir quelles séquences on a déjà vu.
                        Vu.append(False)
                if Vu[target]: #Si on a déjà vu cette sequence, on continue
                    continue
                else:#Sinon on garde la ligne
                    o.write(line)
                    Vu[target] = True
elif len(Arg) == 5 and Arg[3] == "-t": #Cas TE:
    with open(Arg[1],'r') as f:
        with open(Arg[2],'w') as o:
            Vu = set()
            for line in f:
                if line[0] == "#":
                    continue
                target = line.split("\t")[int(Arg[4])] #On récupère la target
                if ';' in target: #On fait un découpage plus fin de la target afin de garder uniquement le gène et non le transcrit
                    target = target.split(';')[0]
                elif ' ' in target:
                    target = target.split(' ')[0]
                if target in Vu: #Si on a déjà vu cette sequence, on continue
                    continue
                else: #Sinon on garde la ligne
                    o.write(line)
                    Vu.add(target[:])
else: #Supprime les lignes identiques se suivant
    with open(Arg[1],'r') as f:
        with open(Arg[2],'w') as o:
            oldline = ""
            for line in f:
                if oldline == "":
                    oldline = line
                    o.write(oldline)
                elif line == oldline:
                    continue
                else:
                    o.write(line)
                    oldline=line