#The goal of this function is to analyse the comp.txt file containing sequences in
# order to compute the consensus sequence of the sequences. Since the component is
# connected, the consensus sequence is defined as the sequences of the path containing
# the heaviest node and 

import sys


Arg = sys.argv[:]

if len(Arg) not in [5]:
    print("Use : " + Arg[0] + "comp.txt graph.edges graph.abundance k")
    exit()

