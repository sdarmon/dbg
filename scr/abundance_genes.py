#The goal of thsi file is to compute the abundance of the genes in the comps
#In input :  -expressed_genes.txt : a file containing the gene names that are expressed
#           -all_unitigs_annotated.txt : a file containing all the unitigs annotated with the gene names

import sys


#Check the input
if len(sys.argv) != 4:
    print("Usage: python3 abundance_genes.py expressed_genes.txt all_unitigs_annotated.txt k")
    sys.exit(1)

k = int(sys.argv[3])

#Reading the expressed genes as a dictionary
expressed_genes = {}
expressed_genes_names = set()
with open(sys.argv[1], 'r') as f:
    for line in f:
        L = line[:-1].split('; ')
        for el in L:
            name = el.split('@')[0]
            if el not in expressed_genes:
                expressed_genes[name] = (0,0) #Size of the gene, number of times it is expressed
            if name not in expressed_genes_names:
                expressed_genes_names.add(el)

#Reading the annotated unitigs on the fly
with open(sys.argv[2], 'r') as f:
    for line in f:
        L = line.split('\t')
        if L[4] != "*":
            genes = L[4].split('; ')
            size = len(genes)
            ab= int(L[8])
            length = len(L[1])

            names_gene = set()
            for gene in genes:
                    gene_name = gene.split('@')[0]
                    names_gene.add(gene_name)
            for gene in genes:
                gene_name = gene.split('@')[0]
                if gene_name in expressed_genes and len(names_gene) == 1:
                    expressed_genes[gene_name] = ( expressed_genes[gene_name][0] + length-k+1, expressed_genes[gene_name][1]+ab)
                    break

#Writing the genes with the abundance into a file_annotated
with open(sys.argv[1].split('.')[0] + "_annotated.txt", 'w') as f_out:
    for gene in expressed_genes_names:
        gene_name = gene.split('@')[0]
        f_out.write(gene + '\t' + str(expressed_genes[gene_name][1]/expressed_genes[gene_name][0]) + '\n')
