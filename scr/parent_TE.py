
import sys

Arg = sys.argv[:]

if len(Arg) not in [7]:
    print("Use : " + Arg[0] + " corrected_nodes.txt Te.fa sorted_count.txt output.txt initial_nodes file.abundance")
    exit()

def rev_comp(seq):
    comp = {'A':'T', 'T':'A', 'C':'G', 'G':'C', 'N':'N'}
    return ''.join([comp[el] for el in seq[::-1]])


#Reading the Te.fa file
Te = {}
name = ''
with open(Arg[2], 'r') as f:
    for line in f:
        if line[0] == '>':
            name = line[1:-1]
            Te[name] = ''
        else:
            Te[name] += line.strip()

#Reading the TE_count.txt file
TE_count = []
with open(Arg[3], 'r') as f:
    for line in f:
        TE_count.append(line.split('\t')[0])

#Define a function that checks if a sequence is a subsequence of a TE
def is_subseq(seq, TE):
    for i in range(len(TE) - len(seq) + 1):
        if TE[i:i+len(seq)] == seq:
            return True
    return False

initial_nodes = {}
#Reading the initial_nodes file
with open(Arg[5], 'r') as f:
    for line in f:
        initial_nodes[line.split('\t')[1]] = line.split('\t')[0]

#Reading the abundance file
abundance = []
with open(Arg[6], 'r') as f:
    for line in f:
        abundance.append(line.split('.')[0])


#Reading the corrected_nodes.txt file and writing the output.txt file
with open(Arg[1], 'r') as f:
    with open(Arg[4], 'w') as f_out:
        for line in f:
            L = line.split('\t')
            seq = L[1]
            rev = rev_comp(seq)
            #looping over the TE the TE_count and find the index of the TE containing the sequence as a subsequence.
            found = False
            found_initial = False
            index_initial = "-1"
            for node,idx in initial_nodes.items():
                if is_subseq(seq, node) or is_subseq(rev, node):
                    found_initial = True
                    index_initial = idx
                    break
            if found_initial:
                ini = "In"
            else:
                ini = "Out"
            for i in range(len(TE_count)):
                if is_subseq(seq, Te[TE_count[i]]) or is_subseq(rev, Te[TE_count[i]]):
                    f_out.write(line[:-1] + '\t' + str(i) + '\t' + TE_count[i] + '\t' + ini + '\t' + abundance[int(index_initial)] + '\n')
                    found = True
                    break
            if not found:
                f_out.write(line[:-1] + '\t'  + "-1" + '\t' + "None" + '\t' + ini + '\t' + abundance[int(index_initial)] + '\n')
