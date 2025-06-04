#The goal is to load the files associated with the comps defined by those rust codes :
#//Saving those values into a file "gene_summary.txt"
#let mut file = std::fs::File::create(format!("{}/gene_summary.txt", output_dir)).unwrap();
#writeln!(file, "{:.2} \t % (Percentage of unitigs covering TE)", percent).unwrap();
#writeln!(file, "{} \t (Number of distinct TE)", number_te).unwrap();
#writeln!(file, "{:.2} \t (Mean distance of the genes)", average_distance_genes).unwrap();
#writeln!(file, "{:.2} \t % (Percentage of exon genes)", proportion_exon_genes).unwrap();
#from analysis_comp import total_length


import matplotlib.pyplot as plt
import numpy as np
import sys
#Check the input
if len(sys.argv) != 3:
    print("Usage: python3 plot_comps.py files_dir nb_comps")
    sys.exit(1)

files_dir = sys.argv[1]
nb_comps = int(sys.argv[2])


#Do some bar graphs to show the distribution of the 4 values types

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
    vu = {}
    for l in range(3,7):
        j = 0
        while j < len(seq) - 2*l :
            s = seq[j:j+l]
            s_rev = rev_comp(s)
            if seq[j+l:j+2*l] == s and seq[j+2*l:j+3*l] == s:
                if j-l >= 0 and s == seq[j-l:j] :
                    if s in vu:
                        vu[s] += l
                    else :
                        vu[s_rev] +=l
                else :
                    if s in vu:
                        vu[s] += 3*l
                    elif s_rev in vu :
                        vu[s_rev] +=3*l
                    else :
                        vu[s] = 3*l
            j += 1
    #Finding the highest sum of vu[*][i] for i in range(6,(encode("TTTTTT")+1))
    m = 0
    seq_m = ''
    for s in vu:
        nb_copy = vu[s]
        if nb_copy > m:
            m = nb_copy
            seq_m = s
    return (m, seq_m, len(seq_m)*m/len(seq))


#Reading the files
percent = []
number_te = []
average_distance_genes = []
proportion_exon_genes = []
polyA = []
microsat = []
for i in range(nb_comps):
    seqs_comp = []
    total_poly = 0
    total_length = 0
    with open("{}/genes_of_comp{}/gene_summary.txt".format(files_dir, i), 'r') as f:
        percent.append(float(f.readline().split('\t')[0]))
        number_te.append(int(f.readline().split('\t')[0]))
        average_distance_genes.append(float(f.readline().split('\t')[0]))
        proportion_exon_genes.append(float(f.readline().split('\t')[0]))

    with open("{}/comp{}_annotated.nodes".format(files_dir, i), 'r') as f:
        for line in f.readlines() :
            seqs_comp.append(line.split('\t')[1])
            total_poly += count_poly(seqs_comp[-1])
            total_length += len(seqs_comp[-1])

    dic_microsat = {}
    for seq in seqs_comp:
        m, seq_m, r = count_microsat(seq)
        if r > 0.5:
            if seq_m in dic_microsat:
                dic_microsat[seq_m] += 1
            elif rev_comp(seq_m) in dic_microsat:
                dic_microsat[rev_comp(seq_m)] += 1
            else:
                dic_microsat[seq_m] = 1
    s = ""
    m_max = -1
    for seq_m in dic_microsat:
        if dic_microsat[seq_m] > m_max:
            m_max = dic_microsat[seq_m]
            s = seq_m

    ratio = m_max / len(seqs_comp)
    polyA_ratio = total_poly / len(seqs_comp)
    if polyA_ratio > 0.5:
        polyA.append(True)
    else:
        polyA.append(False)
    if ratio > 0.20 :
        microsat.append(True)
        print(i, s,m_max, ratio,len(seqs_comp))
    else:
        microsat.append(False)
#Plotting the values

fig, axs = plt.subplots(2, 2)
fig.suptitle('Distribution of the values for the {} first comps'.format(nb_comps))

#Let define four categories, one color for each category
colors = ['green', 'blue', 'red', 'yellow'] # [PolyA,microsat] : False,False = green, False,True = blue, True,False = red, True,True = yellow
#Define a function that map i to the color of the category
def color(i,PolyA,microsat):
    if PolyA[i] and microsat[i]:
        return colors[3]
    elif PolyA[i]:
        return colors[2]
    elif microsat[i]:
        return colors[1]
    else:
        return colors[0]


def label_color(color):
    if color == colors[0]:
        return 'None'
    elif color == colors[1]:
        return 'Microsat'
    elif color == colors[2]:
        return 'A/T Stretches'
    else:
        return 'A/T Stretches and microsat'

# Filter out NaN values
percent2 = [(percent[i],color(i,polyA,microsat)) for i in range(len(percent)) if not np.isnan(percent[i])]
number_te2 = [(number_te[i],color(i,polyA,microsat)) for i in range(len(number_te)) if not np.isnan(number_te[i])]
average_distance_genes2 = [(average_distance_genes[i],color(i,polyA,microsat)) for i in range(len(average_distance_genes)) if not np.isnan(average_distance_genes[i])]
proportion_exon_genes2 = [(proportion_exon_genes[i],color(i,polyA,microsat)) for i in range(len(proportion_exon_genes)) if not np.isnan(proportion_exon_genes[i])]


# Calculate heights for each histogram and convert to percentage, plot the histograms for
# each color category one above the other from green to yellow
percent_heights, percent_bins = np.histogram([x[0] for x in percent2], bins=20)
percent_heights = (percent_heights / nb_comps) * 100

# Initialize a dictionary to store the heights for each color
color_heights = {color: [0] * 20 for color in colors}

# Compute the heights for each color
for i in range(len(percent2)):
    bin_index = 0
    while bin_index + 1 < len(percent_bins) and percent2[i][0] > percent_bins[bin_index + 1]:
        bin_index += 1
    color_heights[percent2[i][1]][bin_index] += 1/nb_comps * 100

# Plot the bars for each color, stacking them
bottom = [0] * len(percent_bins[:-1])
for color in colors:
    axs[0, 0].bar(percent_bins[:-1], color_heights[color], width=np.diff(percent_bins), color=color, align="edge", label=label_color(color), bottom=bottom)
    bottom = [bottom[i] + color_heights[color][i] for i in range(len(bottom))]

axs[0, 0].set_title('Percentage of unitigs covering TE')
axs[0, 0].legend()


#NumberTe case:
b = max(number_te) - min(number_te) + 1
number_te_heights, number_te_bins = np.histogram([x[0] for x in number_te2], bins=b)
number_te_heights = (number_te_heights / nb_comps) * 100

# Initialize a dictionary to store the heights for each color
color_heights = {color: [0] * b for color in colors}

# Compute the heights for each color
for i in range(len(number_te2)):
    bin_index = 0
    while bin_index + 1 < len(number_te_bins) and number_te2[i][0] > number_te_bins[bin_index + 1]:
        bin_index += 1
    color_heights[number_te2[i][1]][bin_index] += 1/nb_comps * 100

# Plot the bars for each color, stacking them

bottom = [0] * len(number_te_bins[:-1])
for color in colors:
    axs[0, 1].bar(number_te_bins[:-1], color_heights[color], width=np.diff(number_te_bins), color=color, align="edge", label=label_color(color), bottom=bottom)
    bottom = [bottom[i] + color_heights[color][i] for i in range(len(bottom))]
axs[0, 1].set_title('Number of distinct TE')
axs[0, 1].legend()

#Average distance genes case:
average_distance_genes_heights, average_distance_genes_bins = np.histogram([x[0] for x in average_distance_genes2], bins=20)
average_distance_genes_heights = (average_distance_genes_heights / nb_comps) * 100

# Initialize a dictionary to store the heights for each color
color_heights = {color: [0] * 20 for color in colors}

# Compute the heights for each color
for i in range(len(average_distance_genes2)):
    bin_index = 0
    while bin_index + 1 < len(average_distance_genes_bins) and average_distance_genes2[i][0] > average_distance_genes_bins[bin_index + 1]:
        bin_index += 1
    color_heights[average_distance_genes2[i][1]][bin_index] += 1/nb_comps * 100

# Plot the bars for each color, stacking them
bottom = [0] * len(average_distance_genes_bins[:-1])
for color in colors:
    axs[1, 0].bar(average_distance_genes_bins[:-1], color_heights[color], width=np.diff(average_distance_genes_bins), color=color, align="edge", label=label_color(color), bottom=bottom)
    bottom = [bottom[i] + color_heights[color][i] for i in range(len(bottom))]
axs[1, 0].set_title('Mean distance of the genes')
axs[1, 0].legend()

#Proportion exon genes case:
proportion_exon_genes_heights, proportion_exon_genes_bins = np.histogram([x[0] for x in proportion_exon_genes2], bins=20)
proportion_exon_genes_heights = (proportion_exon_genes_heights / nb_comps) * 100

# Initialize a dictionary to store the heights for each color
color_heights = {color: [0] * 20 for color in colors}

# Compute the heights for each color
for i in range(len(proportion_exon_genes2)):
    bin_index = 0
    while bin_index + 1 < len(proportion_exon_genes_bins) and proportion_exon_genes2[i][0] > proportion_exon_genes_bins[bin_index + 1]:
        bin_index += 1
    color_heights[proportion_exon_genes2[i][1]][bin_index] += 1/nb_comps * 100

# Plot the bars for each color, stacking them
bottom = [0] * len(proportion_exon_genes_bins[:-1])
for color in colors:
    axs[1, 1].bar(proportion_exon_genes_bins[:-1], color_heights[color], width=np.diff(proportion_exon_genes_bins), color=color, align="edge", label=label_color(color), bottom=bottom)
    bottom = [bottom[i] + color_heights[color][i] for i in range(len(bottom))]

axs[1, 1].set_title('Percentage of exon genes')
axs[1, 1].legend()



# Set labels
for ax in axs.flat: 
    ax.set(xlabel='Values', ylabel='Frequency (%)')

plt.show()
