import math as m

import numpy as np
import matplotlib.pyplot as plt

from numpy import loadtxt
import os
import random

os.chdir("/Users/veronicachang/Desktop/pyThesis")

RM = 15
scale = 3
parNum = 5000000
sig = 10 ** (-5)


def den(r, s):
    return np.exp(-r / s)


def pDF(r, s, Rm):
    up = np.pi * r * den(r, s)
    down = np.pi * s * (s - den(Rm, s) * (Rm + s))
    return up / down


def cDF(r, s, Rm):
    up = np.pi * s * (s - den(r, s) * (r + s))
    down = np.pi * s * (s - den(Rm, s) * (Rm + s))
    return up / down


def invHelper(r, rand, s, Rm):
    cPart = np.pi * s * (s - den(Rm, s) * (Rm + s))
    # left = np.log(s - rand / (m.pi * s * cPart))
    left = np.log(s - rand * cPart / (m.pi * s))
    right = (-r / s) + np.log(r + s)

    return right - left


def invCDF(rand, s, Rm):
    rmin = 0
    rmax = RM
    while (rmax - rmin) / 2 > sig:
        rmid = (rmax + rmin) / 2
        if invHelper(rmid, rand, s, Rm) * invHelper(rmax, rand, s, Rm) < 0:
            rmin = rmid
        else:
            rmax = rmid
    return rmid


def getRsamples(parNum, FILEPATH):
    random_samples = np.array([random.random() for i in range(parNum)])
    samples = [invCDF(rand, scale, RM) for rand in random_samples]
    theta_list = np.array([2 * m.pi * random.random() for i in range(parNum)])

    xis = samples * np.cos(theta_list)
    yis = samples * np.sin(theta_list)

    mass_r = open(FILEPATH, "w")

    for i in range(parNum):
        mass_r.write("%f\t %f\n" % (xis[i], yis[i]))

    mass_r.close()
    print("generation finished")


def plotPdfCdf():
    random_samples = np.array([random.random() for i in range(5000)])
    samples = [invCDF(rand, scale, RM) for rand in random_samples]

    (figure, axes) = plt.subplots()
    xs = np.arange(0., RM, 0.01)
    yes = den(xs, scale)
    yps = pDF(xs, scale, RM)
    yns = cDF(xs, scale, RM)

    ys = [invHelper(i, 0.4, 3, 15) for i in xs]

    plt.plot(xs, yes, label="Column density", color="#C70039", linewidth=2.5)
    plt.plot(xs, yps, label="PDF", color="#7D3C98", linewidth=2.5)
    plt.plot(xs, yns, label="CDF ", color="#3498DB", linewidth=2.5)
    plt.plot(xs, ys, label="CDF monotonicity ", linestyle='--', color="#D35400", linewidth=2.5, )
    plt.ylim(10 ** (-4), 5)
    plt.yscale("log")
    axes.yaxis.set_tick_params(direction='in', which='both')
    axes.xaxis.set_tick_params(direction='in', which='both')

    num_bins = 30
    n, bins, patches = plt.hist(samples, num_bins, facecolor='#D7BDE2', density=True, alpha=0.5)

    plt.title("Exponential Initial Column Surface Density")
    plt.xlabel("R [kpc]")
    plt.ylabel("Value")
    plt.legend()


    return plt.savefig('images/'+ 'E_iniColDenPdfCdf.png', dpi=200)
