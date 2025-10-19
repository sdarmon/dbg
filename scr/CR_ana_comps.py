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
both = 0
UTR_CDS = 0
TE = 0
Nothing = 0

# Process the input files
for i in range(total):
    with open(f"{input_folder}{i}.txt", 'r') as f:
        for line in f:
            L = line.split('\t')
            polA += (1 if float(L[0]) >= 1 else 0)
            microsat += (1 if float(L[1]) >= 0.2 else 0)
            UTR_CDS += (1 if int(L[2]) > 2 else 0)
            TE += (1 if int(L[3]) > 0 else 0)
            #Case otherwise
            if int(L[3]) > 0 :
                print(f"TE")
            elif float(L[0]) >= 0.80 and float(L[1]) >= 0.2:
                both += 1
                print(f"Microsatellite and A/T stretch")
            elif float(L[0]) < 0.80 and float(L[1]) < 0.2:
                Nothing+=1
                print(f"Nothing")
            elif float(L[0]) >= 0.80 and float(L[1]) < 0.2:
                print(f"A/T stretch")
            else :
                print(f"Microsatellite")

        # Sum total for the pie chart
sum_tot = polA + microsat + TE + Nothing - both # UTR_CDS is not included in the pie chart


# Custom function to format the percentages
def func(pct, allvals):
    index = int(round(pct/100.*len(allvals)))
    absolute = allvals[index] if index < len(allvals) else 0
    return "{:.1f}%\n({:d})".format(pct, absolute)

# Data for the pie chart
labels = ['Microsatellite','A/T and microsat', 'A/T stretch',   'TE', 'Unknown']
sizes = [microsat-both, both, polA-both, TE, Nothing]
colors = ['#084D4F','#046668','#007E81','#95C11F','#C11F44']


#Remove the empty sets
for i in range(len(sizes)-1, -1, -1):
    if sizes[i] == 0:
        del sizes[i]
        del labels[i]
        del colors[i]

# Plotting the pie chart
plt.figure(figsize=(8, 8))

def make_autopct(values):
    def my_autopct(pct):
        total = sum(values)
        val = int(round(pct*total/100.0))
        return '{:.1f}%\n({:d})'.format(pct, val)
    return my_autopct

wedges, texts, autotexts =  plt.pie(sizes, labels=labels, colors=colors, autopct=make_autopct(sizes), startangle=140)


# Customize the font properties
for text in texts:
    text.set_fontsize(16)  # Set label font size


for text in autotexts:
    text.set_color('white')
    text.set_fontsize(16)

plt.axis('equal')  # Equal aspect ratio ensures that pie is drawn as a circle.

#plt.title('Distribution of the Components', fontsize=20)
plt.show()
