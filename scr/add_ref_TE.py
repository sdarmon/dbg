#This code read a a comp%i.txt file and a seq_intersectionTE%i.txt file
#and edit the comp%i.txt file to add the TE reference in the seq_intersectionTE%i.txt file for every line
# that is present in the seq_intersectionTE%i.txt file


import sys

Arg = sys.argv[:]

if len(Arg) not in [7]:
    print("Use : " + Arg[0] + " comp%i.txt seq_intersectionRef%i.txt seq_intersectionConsensusDfam%i.txt seq_intersectionConsensusRb%i.txt seq_intergenes$i.txt file.abundance")
    exit()

#Reading intersectionRef%i.txt into a list to get the gene intersection with the comps
dic_index2gene = {}
exon_non_codant = {}
potential_intron = {}
"""
gene_biotype = ""
for info in L[14].split(';'):
    if info.startswith(" gene_biotype"):
        gene_biotype = info.split('"')[1]
        break
                
                
elif gene_biotype == "lncRNA" :
typed = "lncRNA"
"""
with open(Arg[2], 'r') as f:
    for line in f:
        L = line.split('\t')
        gene = L[14].split('"')[1]
        index = int(L[3].split('_')[1])
        position = L[8]
        overlap = int(L[15])
        AS = int(L[16])
        if position in ["CDS","start_codon","stop_codon"]:
            typed="CDS"
        elif position in ["five_prime_utr","three_prime_utr"]:
            typed="UTR"
        elif position in ["exon"]:
            #Case of a potential exon that is not in the CDS
            if index not in exon_non_codant:
                exon_non_codant[index] = [(gene,AS,overlap)]
            elif gene not in exon_non_codant[index]:
                exon_non_codant[index].append((gene,AS,overlap))
            continue
        elif position in ["gene", "transcript"]:
            if index not in potential_intron:
                potential_intron[index] = [(gene,AS,overlap)]
            elif gene not in potential_intron[index]:
                potential_intron[index].append((gene,AS,overlap))
            continue
        else: #Should not happen
            continue
        if index not in dic_index2gene:
            dic_index2gene[index] = [(gene,typed,AS,overlap)]
        elif (gene,typed,AS) not in dic_index2gene[index]:
            dic_index2gene[index].append((gene,typed,AS,overlap))
    #Now loop over the potential exon to add the exon to the gene list if it is not already in the list
    for index in exon_non_codant:
        for (gene,AS,overlap) in exon_non_codant[index]:
            if index not in dic_index2gene :
                dic_index2gene[index] = [(gene,"NCE",AS,overlap)]
            elif gene not in [el[0] for el in dic_index2gene[index]] :
                dic_index2gene[index].append((gene,"NCE",AS,overlap))


    #Now loop over the potential intron to add the intron to the gene list if it is not already in the list
    for index in potential_intron:
        for (gene,AS,overlap) in potential_intron[index]:
            if index not in dic_index2gene :
                dic_index2gene[index] = [(gene,"intron",AS,overlap)]
            elif gene not in [el[0] for el in dic_index2gene[index]] :
                dic_index2gene[index].append((gene,"intron",AS,overlap))



#Reading intersectionConsensus%i.txt into a list to get the TE consensus intersection with the comps
dic_index2consDfam = {}

with open(Arg[3], 'r') as f:
    for line in f:
        L = line.split('\t')
        gene = L[2]
        index = int(L[0].split('_')[1])
        AS = int(L[13].split(':')[2])
        if index not in dic_index2consDfam:
            dic_index2consDfam[index] = [(gene,AS)]
        elif gene not in dic_index2consDfam[index]:
            dic_index2consDfam[index].append((gene,AS))


#Reading intersectionConsensusRB%i.txt into a list to get the TE consensus intersection with the comps
dic_index2consRb = {}

with open(Arg[4], 'r') as f:
    for line in f:
        L = line.split('\t')
        if len(L) < 2 :
            continue
        gene = L[2]
        index = int(L[0].split('_')[1])
        if index not in dic_index2consRb:
            dic_index2consRb[index] = [gene]
        elif gene not in dic_index2consRb[index]:
            dic_index2consRb[index].append(gene)


#Reading intersectionGenes%i.txt into a list to get the gene intersection with the comps$i
dic_index2intergene = {}
with open(Arg[5], 'r') as f:
    for line in f:
        L = line.split('\t')
        if len(L) < 2 :
            continue
        index = int(L[0].split('_')[1])
        chromosome = L[2]
        if index not in dic_index2intergene:
            dic_index2intergene[index] = [chromosome]
        elif chromosome not in dic_index2intergene[index]:
            dic_index2intergene[index].append(chromosome)

#Reading the abundance file
abundance = []
with open(Arg[6], 'r') as f:
    for line in f:
        abundance.append(line.split('.')[0])


#Reading every line of comp%i.txt and adding the TE reference if the index of the line is in dic_index2TE
#Writing the new line in a new file comp%i.txt
i = 0
with open(Arg[1], 'r') as f:
    with open(Arg[1].split('.')[0] + "_annotated.nodes", 'w') as f_out:
        for line in f:
            L = line.split('\t')
            node = int(L[0])
            AS=0
            AS_dfam=0
            if i in dic_index2gene:
                str_gene = "; ".join([el[0] + "@" + el[1] + "%" + str(el[3]) for el in dic_index2gene[i]])
                AS = max(el[2] for el in dic_index2gene[i])
            else:
                str_gene = "*"

            if i  in dic_index2consDfam:
                str_consDfam = "; ".join([el[0] for el in dic_index2consDfam[i]])
                AS_dfam = max(el[1] for el in dic_index2consDfam[i])
            else:
                str_consDfam = "*"
            if i  in dic_index2consRb:
                str_consRb = "; ".join(dic_index2consRb[i])
            else:
                str_consRb = "*"
            if i in dic_index2intergene:
                str_intergene = "; ".join(dic_index2intergene[i])
            else:
                str_intergene = "*"
            L=line[:-1].split("\t")
            ind=L[0]
            seq=L[1]
            if len(L) == 3:
                weight= L[2]
                d = "0"
            else:
                d = L[2]
                weight = L[3]
            f_out.write(ind +"\t" + seq + "\t" + d + "\t" + weight + "\t" + str_gene + "\t" + str_consDfam  + "\t" + str_consRb  + "\t" + str_intergene  + "\t" +  abundance[int(L[0])] + "\t" + str(AS) + "\t" + str(AS_dfam) + "\n")

            i += 1
