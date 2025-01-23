#The goal of this file is to compare the weight of the same nodes in two different files
#The files are given as arguments and it should output a graph on which the x-axis is the weight of the nodes of
# the first file and the y-axis is the weight of nodes of the second file

import sys
import matplotlib.pyplot as plt
import numpy as np

#Check the input
if len(sys.argv) != 5:
    print("Usage: python3 dist_comparaison.py file1 file2 comp_prefix nb_comp")
    sys.exit(1)

#Reading the comp files
n = int(sys.argv[4])
comps= [[] for i in range(n)]
for i in range(n):
    with open(sys.argv[3] + str(i) + ".txt", 'r') as f:
        for line in f:
            L = line.split()
            comps[i].append(int(L[0]))


#Read the files the nodes are written on one line as id seq weight
nodesX = []
nodesY = []
with open(sys.argv[1], 'r') as f:
    for line in f:
        L = line.split()
        nodesX.append(int(L[2]))
with open(sys.argv[2], 'r') as f:
    for line in f:
        L = line.split()
        nodesY.append(int(L[2]))

#Compute the correlation coefficient
correlation = np.corrcoef(nodesX, nodesY)[0, 1]
print("Correlation coefficient: " + str(correlation))

#Only keep the nodes that have the most weight (top 1%)
X = []
Y = []
sortedX = sorted(nodesX, reverse=True)
top = sortedX[int(len(sortedX) * 0.01)]
for i in range(len(nodesX)):
    if nodesX[i] >= top:
        X.append(nodesX[i])
        Y.append(nodesY[i])


#plot the graph
plt.scatter(X,Y, color='black')
plt.plot(np.unique(X), np.poly1d(np.polyfit(X, Y, 1))(np.unique(X)), color='red')

#For every nodes inside the comp files, plot them (one color per comp)
for i in range(n):
    Xi = []
    Yi = []
    for j in comps[i]:
        Xi.append(nodesX[j])
        Yi.append(nodesY[j])
    plt.scatter(Xi,Yi, label="comp" + str(i))


plt.xlabel(sys.argv[1])
plt.ylabel(sys.argv[2])
plt.title("Comparaison of the weight of the nodes")
plt.show()