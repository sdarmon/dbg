#The goal of this function is to analyse the comp.txt file containing sequences in order to compute the
#number of poly(A) (represented as five consecutive A's) in the sequences; the ratio of microsatellites
#(repeated sequences of 1 to 6 nucleotides) in the sequences; and the annotated genes that are intersecting
# with the sequences.
import sys
from math import log

Arg = sys.argv[:]

if len(Arg) not in [5,6,7,8]:
    print("Use : " + Arg[0] + " comp.fq intersection_gene.txt intersection_TE.txt output.txt [--o.s. sep] [--short]")
    exit()



#Function that count how many poly(A) a sequence has
def count_poly(seq):
    compt = 0
    for i in range(len(seq) - 4):
        if seq[i:i+5] == 'AAAAA' or seq[i:i+5] == 'TTTTT':
            compt += 1
    return compt

#Function that count the ratio of C/G in a sequence
def count_CG(seq):
    countCG = 0
    countAT = 0
    for el in seq:
        if el == 'C' or el == 'G':
            countCG += 1
        elif el == 'A' or el == 'T':
            countAT += 1
    return countCG / (countCG + countAT)

#Function that encodes a sequence in a binary format (blank=0, A=1, C=2, G=3, T=4)
def encode(seq):
    enc = {'A':1, 'C':2, 'G':3, 'T':4}
    s = 0
    for i in range(len(seq)):
        s = s*5 + enc[seq[i]]
    return s

#Function that decodes a sequence from a binary format to a nucleotide sequence
def decode(s):
    dec = {0:' ', 1:'A', 2:'C', 3:'G', 4:'T'}
    seq = ''
    while s > 0:
        if s%5 == 0:
            return ''
        seq = dec[s%5] + seq
        s = s//5
    return seq

#Function that checks if a sequence is a microsatellite
def is_a_microsat(seq):
    for i in range(1,len(seq)//2+1):
        if len(seq) % i != 0:
            continue
        if seq == seq[:i] * (len(seq) // i):
            return True
    return False

#Function that computes the reverse complement of a sequence
def rev_comp(seq):
    comp = {'A':'T', 'T':'A', 'C':'G', 'G':'C', 'N':'N'}
    return ''.join([comp[el] for el in seq[::-1]])

#Function that counts the number of microsatellites in a sequence
def count_microsat(seq):
    vu = [[0 for _ in range(encode("TTTTT")+1)] for _ in range(len(seq))]
    for i in range(6,(encode("TTTTT")+1)):
        s = decode(i)
        if s == '' or is_a_microsat(s):
            continue
        l = len(s)
        for j in range(len(seq) - l):
            if seq[j:j+l] == s and seq[j+l:j+2*l] == s:
                for k in range(j,j+2*l):
                    vu[k][i] = 1
                    i_rev = encode(rev_comp(s))
                    vu[k][i_rev] = 1
    #Finding the highest sum of vu[*][i] for i in range(6,(encode("TTTTT")+1))
    m = 0
    seq_m = ''
    for i in range(6,(encode("TTTTT")+1)):
        s = decode(i)
        if s == '' or is_a_microsat(s):
            continue
        nb_copy = sum([vu[j][i] for j in range(len(seq))]) / len(s)
        if nb_copy > m:
            m = nb_copy
            seq_m = s

    #Finding the positions covered by a microsatellite
    pos_vu = [0 for _ in range(len(seq))]
    for i in range(len(seq)):
        for j in range(6,(encode("TTTTT")+1)):
            if vu[i][j] == 1:
                pos_vu[i] = 1
    return (m, seq_m, sum(pos_vu)/len(seq))

#Reading the genes from the gene.gtf file
genes = set()
match_prot = set()
UTR = set()
CDS = set()
intron = set()

if len(Arg) >= 7 :
    sep = Arg[6]
else:
    sep = "split('\t')[20]"
func = eval(f"lambda x: x.{sep}")


#Reading the gene.gtf file
with open(Arg[2], 'r') as f:
    for line in f:
        n = len(line.split('\t'))
        if n < 21:
            continue
        target = func(line)
        if "matches_ref_protein \"True\";" in line:
            match_prot.add(target)
        if "UTR" or "five_prime_utr" or "three_prime_utr" in line :
            UTR.add(target)
        if "CDS" or "start_codon" or "stop_codon" in line:
            CDS.add(target)
        genes.add(target)

#Reading the transposable elements from the aligned sequences file
transpo = set()
with open(Arg[3], 'r') as f:
    for line in f:
        target = line.split('\t')[2]
        if target != "*":
            transpo.add(target)

#Function that computes the entropy of the n-mers in the sequences
def entropy(seqs,n):
    count = [0 for _ in range(5**n)]
    for seq in seqs:
        for i in range(len(seq) - n + 1):
            s = encode(seq[i:i+n])
            count[s] += 1
    total = sum(count)
    prob = [el/total for el in count]
    return -sum([el * log(el,2) for el in prob if el > 0])

#Reading the sequences from the comp.txt file and compute the average number of poly(A) and the ratio of microsatellites
seqs = []
total_poly = 0
total_length = 0
with open(Arg[1], 'r') as f:
    for line in f:
        if len(line) < 2:
            break
        if line[0] == '>':
            continue
        seqs.append(line[:-1])
        total_poly += count_poly(line[:-1])
        total_length += len(line[:-1])
m, seq_m, r = count_microsat(' '.join(seqs))

#We define the introns as the sequences that are not in the UTR nor the CDS but still are in a gene
intron = genes - UTR - CDS

#writing the results in the output file
if len(Arg) == 5 or len(Arg) == 7:
    with open(Arg[4], 'w') as f:
        f.write("Average number of poly(A) : " + str(total_poly/len(seqs)) + "\n")
        f.write("Entropy of the 6-mers : " + str(entropy(seqs,6)) + "\n")
        f.write("Ratio of CG : " + str(count_CG(' '.join(seqs))) + "\n")
        f.write("Ratio of microsatellites : " + str(r) + "\n")
        f.write("Most frequent microsatellite : " + seq_m + " with " + str(m) + " copies, covering " + str(m*len(seq_m)/total_length*100) + "% of the sequences\n")
        f.write(str(len(genes)) + " Genes intersecting with the sequences : " + str(genes) + "\n\n")
        f.write(str(len(match_prot)) + " Genes intersecting with the sequences that match a protein : " + str(match_prot) + "\n\n")
        f.write(str(len(UTR)) + " Genes intersecting with the sequences that are in UTR : " + str(UTR) + "\n\n")
        f.write(str(len(CDS)) + " Genes intersecting with the sequences that are in CDS : " + str(CDS) + "\n\n")
        f.write(str(len(intron)) + " Genes intersecting with the sequences that are in intron : " + str(intron) + "\n\n")
        f.write(str(len(transpo)) + " Transposable elements intersecting with the sequences : " + str(transpo) + "\n\n")

else:
    with open(Arg[4], 'w') as f:
        f.write(str(total_poly/len(seqs)) + "\t")
        f.write(str(r) + "\t")
        f.write(str(len(UTR)+len(CDS)) + "\t")
        f.write(str(len(transpo)) + "\n")

