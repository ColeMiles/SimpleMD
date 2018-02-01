#!/usr/bin/python3

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

df = pd.read_csv("../data/summary.dat", sep="\s+")
plt.plot(df["time"], df["temperature"], label=r"$\langle E_K \rangle$")

plt.plot((0.0, df.tail(1)["time"]), (1.0, 1.0), "r--", label="Ideal temp")

plt.xlabel("$t$")
plt.ylabel(r"$T = \frac{2}{3}\langle E_K \rangle$")
plt.title("Temp vs. Time")
plt.legend()
plt.show()