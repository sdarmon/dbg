
import sys

Arg = sys.argv[:]

if len(Arg) not in [4]:
    print("Use : " + Arg[0] + " allNodes.txt Te.fa TeNodes_corrected.txt")
    exit()


#Function that returns the reverse complement of a sequence
def rev_comp(seq):
    comp = {'A':'T', 'T':'A', 'C':'G', 'G':'C', 'N':'N'}
    return ''.join([comp[el] for el in seq[::-1]])


#Reading the TE.fa file and storing the 41-mers of the TEs
Te = {}
name = ''
with open(Arg[2], 'r') as f:
    for line in f:
        if line[0] == '>':
            name = line[1:-1]
            Te[name] = ''
        else:
            Te[name] += line.strip()


#Define a function that checks if a sequence is a subsequence of a TE
def is_subseq(seq, TE):
    for i in range(len(TE) - len(seq) + 1):
        if TE[i:i+len(seq)] == seq:
            return True
    return False


#Reading the allNodes.txt file and writing the corrected file
with open(Arg[1], 'r') as f:
    with open(Arg[3], 'w') as f_out:
        for line in f:
            L = line.split('\t')
            #Checking if the L[1] is in the TEnodes
            for TE in Te:
                if is_subseq(L[1], Te[TE]) or is_subseq(rev_comp(L[1]), Te[TE]):
                    f_out.write(line)
                    break