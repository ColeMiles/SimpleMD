#!/usr/bin/python3

import subprocess
import os
import itertools

if os.getcwd().endswith("scripts"):
    os.chdir("..")

p_nums = [125, 512, 1000, 1728]
temps = [0.5, 0.75, 1.0, 1.25, 1.5, 1.75, 2.0]
volume = 9.0 ** 3
steps = 20000

cmd = "./driver"
for (pn, temp) in itertools.product(p_nums, temps):
    print("Running sim with temp = {} and {} particles for {} steps."
            .format(temp, pn, steps))
    density = pn / volume

    ghistfile = "gtests/T{}_p{:.2}.dat".format(temp, density) 
    cmdfull = [cmd, "-n", "20000", "-g", ghistfile, "-T", str(temp), "-N", str(pn), "-l"]
    subprocess.run(cmdfull, stdout=subprocess.DEVNULL, check=True)