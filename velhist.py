import matplotlib.pyplot as plt
import numpy as np
from ast import literal_eval

with open("boltzmann.dat", "r") as file:
    #Skip first line
    file.readline()
    for line in file:
        if line == "\n":
            continue
        l_split = line.strip().split("\t")
        time = l_split[0]
        bins, heights = [], []
        for hist_bin in l_split[1:]:
            tuple = literal_eval(hist_bin)
            bins.append(tuple[0])
            heights.append(tuple[1])
        plt.bar(bins, heights, width=bins[1] - bins[0], align="edge", 
                               alpha=0.5, label="$t={}$".format(time))
    plt.legend()
    plt.show()