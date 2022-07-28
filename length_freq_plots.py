import os
import csv
import matplotlib.pyplot as plt
import matplotlib as mpl

mpl.use('macosx')

#chosee to which chromosome you would like to do your experiments
chr=1
pet=3
technology="chiapet"
cell_line = "GM19238_merged_replicates"

path_to_locdir = f"/Users/aleksander.ostruk/Documents/asia_master/"
filename="len_freq_%s_chr%s_pet%s.csv" %(cell_line, chr, pet)
filepath=os.path.join(os.getcwd(), path_to_locdir, technology, cell_line, "length_analysis", filename)

length = []
frequency = []

with open(filepath, 'r') as csvfile:
    rows = list(csv.reader(csvfile, delimiter=' '))      # Read all rows into a list
    for row in rows:
        length.append(int(row[0]))
        frequency.append(int(row[1]))

x = length
y = frequency

def get_name(technology):
    if technology=="hichip":
        prop_name="HiChIP"
    elif technology=="chiapet" :
        prop_name="ChiA-PET"
    return(prop_name)


def generate_scatter_plot(x, y):
    plt.scatter(x, y, color='g', label=f"chr{chr}_pet{pet}", s=5)
    plt.title(f"Length analysis for {cell_line}, {get_name(technology)}")
    plt.autoscale(enable=True, axis='both', tight=None)
    plt.xscale('log')
    plt.xlabel('length [log bp]')
    plt.ylabel('Numbers of PETs')
    plt.legend()
    plt.show()

generate_scatter_plot(x,y)






