#!/usr/bin/python3

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

df = pd.read_csv("../data/summary.dat", sep="\s+")
plt.plot(df["time"], df["avg_e_kin"], label=r"$\langle E_K \rangle$")
plt.plot(df["time"], df["avg_e_pot"], label=r"$\langle E_U \rangle$")
plt.plot(df["time"], df["avg_e_tot"], label=r"$\langle E_{tot} \rangle$")

plt.xlabel("$t$")
plt.ylabel(r"$\langle E \rangle$")
plt.title("Energies vs. Time")
plt.legend()
plt.show()
