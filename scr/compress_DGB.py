# The goal of this script is to compress the new De Bruijn Graph.
# Some edges have been removed so we need to compress the paths of the graph, that is to say,
# the consecutive unitigs of degree 2.
# The compression is done by compacting those unitigs.
# The script takes 3 arguments:
# - the nodes of the graph : each line is the inter index of the node and the sequence of the node
# - the edges of the graph : each line is the index of the source node, the index of the target node and
# direction of that edges (RR, RF, FR or FF)
# - the output file

import sys

#define a function that do a sorted insertion in a list by dichotomy
def insert_sorted(L, e, label):
    if len(L) == 0:
        L.append((e,label))
        return
    if e < L[0]:
        L.insert(0, (e,label))
        return
    if e > L[-1]:
        L.append((e,label))
        return
    a = 0
    b = len(L) - 1
    while b - a > 1:
        m = (a + b) // 2
        if e == L[m]:
            return
        if e < L[m]:
            b = m
        else:
            a = m
    L.insert(b, (e,label))
    return

# Define a function that a sorted deletion in a list by dichotomy
def delete_sorted(L, e):
    if len(L) == 0:
        return
    if e < L[0][0]:
        return
    if e > L[-1][0]:
        return
    a = 0
    b = len(L) - 1
    while b - a > 1:
        m = (a + b) // 2
        if e == L[m][0]:
            L.pop(m)
            return
        if e < L[m][0]:
            b = m
        else:
            a = m
    if e == L[a][0]:
        L.pop(a)
    if e == L[b][0]:
        L.pop(b)
    return

Arg = sys.argv[:]
if len(Arg) not in [4]:
    print("Use : " + Arg[0] + " nodes edges output")
    exit()

# Define V as the list of nodes
V = []
# Define E as the list of edges
E = []

# Read the nodes
with open(Arg[1], 'r') as f:
    for line in f:
        V.append(line.split()[1])

# Read the edges
with open(Arg[2], 'r') as f:
    for line in f:
        E.append(line.split())

# Define the sorted adjacency list of the directed egdes

AdjE = [[] for i in range(len(V))]
for e in E:
    insert_sorted(AdjE[int(e[0])], int(e[1]), e[2])

# Now we scan every node of the graph. If the node has a degree of , we compress the path

