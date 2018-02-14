#!/usr/bin/python3

import matplotlib.pyplot as plt
import numpy as np
import os
import re
from ast import literal_eval

def plot(filename):
    fig = plt.figure()
    print(filename)
    match = re.match("T(.*)_p(.*)\.dat", filename)
    temp = float(match.group(1))
    rho = float(match.group(2))

    with open("data/gtests/" + filename, "r") as file:
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
            plt.plot(bins, heights)        

            # plt.bar(bins, heights, width=bins[1] - bins[0], align="edge", 
            #                        alpha=0.5, label="Data")

        plt.xlabel("$r$")
        plt.ylabel("$g$")
        plt.title(r"$g$ Distribution, $T = {}, \rho = {}$".format(temp, rho))
        fig.savefig("data/gtestimgs/" + "T{}_p{}.png".format(temp, rho))
    plt.close()

if __name__ == "__main__":
    if os.getcwd().endswith("plotting"):
        os.chdir("..")
    files = os.listdir("data/gtests")
    for file in files:
        plot(file)
