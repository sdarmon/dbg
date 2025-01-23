#The goal of this function is to read a gff file, finding the genes and the exons annotations and then
#to add the introns to the gff file as new annotations.

import sys

Arg = sys.argv[:]

if len(Arg) not in [3]:
    print("Use : " + Arg[0] + " input.gff output.gff")
    exit()

genes = [] #List of genes. Each gene is a list composed of the line name then tuples of start and end of the exons.

with open(Arg[1], 'r') as f:
    for line in f:
        if len(line) < 2:
            break
        L = line.split('\t')
        if L[2] == 'gene':
            genes.append([line])
        elif L[2] == 'exon':
            genes[-1].append((int(L[3]),int(L[4])))


#Let's compute the introns
introns = [] #List of introns. Each intron is a list containing the line name then the start and end of the intron.

for gene in genes:
    #Sorting the exons
    exons = sorted(gene[1:], key=lambda x: x[0])
    for i in range(1,len(exons)):
        introns.append([gene[0],exons[i-1][1],exons[i][0]])

#Writing the output file while reformating the introns names
with open(Arg[2], 'w') as f:
    for gene in genes:
        for el in gene:
            f.write(el[0])
    for intron in introns:
        f.write(intron[0].split('\t')[0] + '\t' + intron[0].split('\t')[1] + '\tintron\t' + str(intron[1]) + '\t' + str(intron[2]) + '\t' + '\t'.join(intron[0].split('\t')[5:]))

