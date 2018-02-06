#!/usr/bin/python3

import matplotlib.pyplot as plt
import numpy as np
from ast import literal_eval

with open("../data/ghist.dat", "r") as file:
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
                               alpha=0.5, label="Data")

    plt.xlabel("$g$")
    plt.ylabel("Rel. Frequency")
    plt.title(r"$g$ Distribution, Langevin, $\rho = 0.75$")
    plt.legend()
    plt.show()