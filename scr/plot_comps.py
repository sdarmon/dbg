#The goal is to load the files associated with the comps defined by those rust codes :
#//Saving those values into a file "gene_summary.txt"
#let mut file = std::fs::File::create(format!("{}/gene_summary.txt", output_dir)).unwrap();
#writeln!(file, "{:.2} \t % (Percentage of unitigs covering TE)", percent).unwrap();
#writeln!(file, "{} \t (Number of distinct TE)", number_te).unwrap();
#writeln!(file, "{:.2} \t (Mean distance of the genes)", average_distance_genes).unwrap();
#writeln!(file, "{:.2} \t % (Percentage of exon genes)", proportion_exon_genes).unwrap();

#Do some bar graphs to show the distribution of the 4 values types


import sys
#Check the input
if len(sys.argv) != 3:
    print("Usage: python3 plot_comps.py files_dir nb_comps")
    sys.exit(1)

files_dir = sys.argv[1]
nb_comps = int(sys.argv[2])

#Reading the files
percent = []
number_te = []
average_distance_genes = []
proportion_exon_genes = []

for i in range(nb_comps):
    with open("{}/genes_of_comp{}/gene_summary.txt".format(files_dir, i), 'r') as f:
        percent.append(float(f.readline().split('\t')[0]))
        number_te.append(int(f.readline().split('\t')[0]))
        average_distance_genes.append(float(f.readline().split('\t')[0]))
        proportion_exon_genes.append(float(f.readline().split('\t')[0]))

#Plotting the values
import matplotlib.pyplot as plt
import numpy as np

fig, axs = plt.subplots(2, 2)
fig.suptitle('Distribution of the values for the {} comps'.format(nb_comps))

# Filter out NaN values
percent = [x for x in percent if not np.isnan(x)]
number_te = [x for x in number_te if not np.isnan(x)]
average_distance_genes = [x for x in average_distance_genes if not np.isnan(x)]
proportion_exon_genes = [x for x in proportion_exon_genes if not np.isnan(x)]

# Calculate heights for each histogram and convert to percentage
percent_heights, percent_bins = np.histogram(percent, bins=20)
percent_heights = (percent_heights / nb_comps) * 100
axs[0, 0].bar(percent_bins[:-1], percent_heights, width=np.diff(percent_bins), color="green", align="edge")
axs[0, 0].set_title('Percentage of unitigs covering TE')

number_te_heights, number_te_bins = np.histogram(number_te, bins=range(min(number_te), max(number_te) + 1, 1))
number_te_heights = (number_te_heights / nb_comps) * 100
axs[0, 1].bar(number_te_bins[:-1], number_te_heights, width=np.diff(number_te_bins), color="green", align="edge")
axs[0, 1].set_title('Number of distinct TE')

average_distance_genes_heights, average_distance_genes_bins = np.histogram(average_distance_genes, bins=20)
average_distance_genes_heights = (average_distance_genes_heights / nb_comps) * 100
axs[1, 0].bar(average_distance_genes_bins[:-1], average_distance_genes_heights, width=np.diff(average_distance_genes_bins), color="green", align="edge")
axs[1, 0].set_title('Mean distance of the genes')

proportion_exon_genes_heights, proportion_exon_genes_bins = np.histogram(proportion_exon_genes, bins=20)
proportion_exon_genes_heights = (proportion_exon_genes_heights / nb_comps) * 100
axs[1, 1].bar(proportion_exon_genes_bins[:-1], proportion_exon_genes_heights, width=np.diff(proportion_exon_genes_bins), color="green", align="edge")
axs[1, 1].set_title('Percentage of exon genes')

# Set labels
for ax in axs.flat: 
    ax.set(xlabel='Values', ylabel='Frequency (%)')

plt.show()
