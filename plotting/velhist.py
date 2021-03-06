#!/usr/bin/python3

import matplotlib.pyplot as plt
import numpy as np
from ast import literal_eval

temp = 1.0

def maxwell(v):
    return np.power(2 * np.pi * temp, -3.0 / 2.0) * 4 * np.pi * v * v * np.exp(-v * v / (2.0 * temp))

with open("../data/boltzmann.dat", "r") as file:
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

    #Plot the Maxwell-Boltzman distribution from [0.0, 2.0]
    vspace = np.linspace(0.0, 5.0, num=50)
    plt.xlabel("$v$")
    plt.ylabel("Rel. Frequency")
    plt.plot(vspace, maxwell(vspace), "k-", label="Maxwell-Boltzman")
    plt.title(r"Velocity Distribution, Langevin, $\rho = 0.75$")
    plt.legend()
    plt.show()
