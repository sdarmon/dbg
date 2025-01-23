import sys
import pandas as pd
import matplotlib.pyplot as plt

# Argument processing
Arg = sys.argv[:]
if len(Arg) != 3:
    print("Use : " + Arg[0] + " input_folder total")
    exit()

input_folder = Arg[1]
total = int(Arg[2])

# Initialize counts
polA = 0
microsat = 0
UTR_CDS = 0
TE = 0

# Process the input files
for i in range(total):
    with open(f"{input_folder}{i}.txt", 'r') as f:
        for line in f:
            L = line.split('\t')
            polA += (1 if float(L[0]) >= 1 else 0)
            microsat += (1 if float(L[1]) >= 0.2 else 0)
            UTR_CDS += (1 if int(L[2]) > 2 else 0)
            TE += (1 if int(L[3]) > 0 else 0)

# Sum total for the pie chart
sum_tot = polA + microsat + TE # UTR_CDS is not included in the pie chart


# Custom function to format the percentages
def func(pct, allvals):
    absolute = int(pct/100.*sum(allvals))
    return "{:.1f}%\n({:d})".format(pct, absolute)


# Data for the pie chart
labels = ['A/T stretch', 'Micro-Sat',  'TE']
sizes = [polA, microsat, TE]
colors = ['#007E81','#084D4F','#95C11F']

# Plotting the pie chart
plt.figure(figsize=(8, 8))
wedges, texts, autotexts =  plt.pie(sizes, labels=labels, colors=colors, autopct=lambda pct: func(pct, sizes), startangle=140)


# Customize the font properties
for text in texts:
    text.set_fontsize(16)  # Set label font size


for text in autotexts:
    text.set_color('white')
    text.set_fontsize(16)

plt.axis('equal')  # Equal aspect ratio ensures that pie is drawn as a circle.

plt.title('Distribution of the Components', fontsize=20)
plt.show()
