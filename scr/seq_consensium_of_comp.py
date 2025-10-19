#The goal of this function is to analyse the comp.txt file containing sequences in
# order to compute the consensus sequence of the sequences. Since the component is
# connected, the consensus sequence is defined as the sequences of the path containing
# the heaviest node and 

import sys


Arg = sys.argv[:]

if len(Arg) not in [4]:
    print("Use : " + Arg[0] + " comp_annotated.nodes graph.edges k")
    exit()


#Initialization of the dictionaries
comp = {}
max_node = -1
max_weight = -1
#Reading the comp_annotated.nodes file into a list
#Line format : ind +"\t" + seq + "\t" + str(d) + "\t" + weight + "\t" + str_gene + "\t" + str_consDfam  + "\t" + str_consRb  + "\t" + str_intergene  + "\t" +  abundance[int(L[0])] + "\t" + str(AS) + "\t" + str(AS_dfam) + "\n"
with open(Arg[1], 'r') as f:
    for line in f:
        L=line[:-1].split('\t')
        ind = int(L[0])
        seq = L[1]
        w = int(L[8])
        ab = int(L[8])
        comp[ind] = [seq, w, ab]
        if w > max_weight:
            max_weight = w
            max_node = ind

#Reading the graph.edges file to get the neighbours of each node
#edge format : node1 + "\t" + node2 + "\t" + type + "\n"
graph = {}
with open(Arg[2], 'r') as f:
    for line in f:
        L = line[:-1].split('\t')
        node1 = int(L[0])
        node2 = int(L[1])
        type = L[2][0:2]  # Remove the newline character at the end

        if node1 not in graph:
            graph[node1] = []

        #Checking if the node1 and 2 are in the comp keys
        if node1 in comp and node2 in comp:
            graph[node1].append((node2,type))

        #The reverse edge should always be present in the edge file


nodes_seen = set()
nodes_seen.add(max_node)
forward_path = []
way = 'F'
head = max_node
# Start from the node with the maximum weight

while True:
    if head not in graph:
        break
    neighbours = graph[head]
    next_node = -1
    best_ab = -1
    new_way = 'N'

    for neighbour, edge_type in neighbours:
        if neighbour not in nodes_seen and edge_type[0] == way:
            if comp[neighbour][2] > best_ab:
                next_node = neighbour
                best_ab = comp[neighbour][2]
                new_way = edge_type[1]
    if next_node == -1:
        break

    forward_path.append((next_node, new_way))
    head = next_node
    nodes_seen.add(next_node)

#Idem for the reverse path
reverse_path = []
way = 'R'
head = max_node

while True:
    if head not in graph:
        break
    neighbours = graph[head]
    next_node = -1
    best_ab = -1
    new_way = 'N'

    for neighbour, edge_type in neighbours:
        if neighbour not in nodes_seen and edge_type[0] == way:
            if comp[neighbour][2] > best_ab:
                next_node = neighbour
                best_ab = comp[neighbour][2]
                new_way = edge_type[1]
    if next_node == -1:
        break

    reverse_path.append((next_node, new_way))
    head = next_node
    nodes_seen.add(next_node)

#Now we have the forward and reverse paths, we can compute the consensus sequence

def reverse_complement(seq):
    complement = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
    return ''.join(complement[base] for base in seq[::-1])

def cat(seq1, seq2, way2, k):
    #We suppose that seq1 has a forward orientation
    if way2 == 'F':
        s2 = seq2
    else:
        s2 = reverse_complement(seq2)

    return seq1 + s2[k-1:]

#Compute the consensus sequence
consensus_seq = comp[max_node][0]
k = int(Arg[3])
for node, way in forward_path:
    consensus_seq = cat(consensus_seq, comp[node][0],way, k)

consensus_seq= reverse_complement(consensus_seq)

for node, way in reverse_path:
    consensus_seq = cat(consensus_seq, comp[node][0],way, k)

#Output the consensus sequence
print(consensus_seq)
print("Stats: Total nb nodes :" + str(len(forward_path)+ len(reverse_path) + 1))