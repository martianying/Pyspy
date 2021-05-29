# ----------------------------------------------------------------------------#
# MODULE: iniDensity                                                          #
#                                                                             #
# CONTENT:                                                                    #
#      generate initial gas particles density distribution to follow an       #
#   exponential distribution. Proud to had coded this in PHANTOM in FORTRAN   #
#   language!                                                                 #
#                                                                             #
# AIM:                                                                        #
#      1) generate samples of desired distribution                            #
#      2) plot the scatter plot for samples together with the theoretical     #
#      curve to check its validity.                                           #
#                                                                             #
# DEPENDENCIES: phycon                                                        #
# ----------------------------------------------------------------------------#
import os
import random
import matplotlib.pyplot as plt
from src.phycon import *

SCALE = 3
parNum = 5000000
SIG = 10 ** (-5)

def den(r, s):
    """
    density function
    """
    return np.exp(-r / s)


def pDF(r, s, Rm):
    """
    probability density function
    """
    up = np.pi * r * den(r, s)
    down = np.pi * s * (s - den(Rm, s) * (Rm + s))
    return up / down


def cDF(r, s, Rm):
    """
    cumulative density function
    """
    up = np.pi * s * (s - den(r, s) * (r + s))
    down = np.pi * s * (s - den(Rm, s) * (Rm + s))
    return up / down


def invHelper(r, rand, s, Rm):
    """
    Helper function to calculate inverse cumulative density function
    """
    cPart = np.pi * s * (s - den(Rm, s) * (Rm + s))
    left = np.log(s - rand * cPart / (m.pi * s))
    right = (-r / s) + np.log(r + s)

    return right - left


def invCDF(rand, s, Rm):
    """
    inverse cumulative function
    """
    rmin = 0
    rmax = RM
    while (rmax - rmin) / 2 > SIG:
        rmid = (rmax + rmin) / 2
        if invHelper(rmid, rand, s, Rm) * invHelper(rmax, rand, s, Rm) < 0:
            rmin = rmid
        else:
            rmax = rmid
    return rmid


def generateAndSaveRSamples(parNum, FILEPATH):
    """
    HElPER function
    generate and save samples in a new file

    Arguments:
    parNum:   number of sampling particles
    FILEPATH: set the path for the file of samples

    Returns:
    file full of samples at location FILEPATH
    """
    random_samples = np.array([random.random() for i in range(parNum)])
    samples = [invCDF(rand, SCALE, RM) for rand in random_samples]
    theta_list = np.array([2 * m.pi * random.random() for i in range(parNum)])

    # calculate the postions
    xis = samples * np.cos(theta_list)
    yis = samples * np.sin(theta_list)

    mass_r = open(FILEPATH, "w")

    for i in range(parNum):
        mass_r.write("%f\t %f\n" % (xis[i], yis[i]))
    mass_r.close()

    print("generation finished")

def getAndPlotRsamples(parNum):
    """
    same as the function name
    """
    plotName = 'Exp Setup Python Scatter Plot.png'
    random_samples = np.array([random.random() for i in range(parNum)])
    samples = [invCDF(rand, SCALE, RM) for rand in random_samples]
    theta_list = np.array([2 * m.pi * random.random() for i in range(parNum)])

    xis = samples * np.cos(theta_list)
    yis = samples * np.sin(theta_list)

    fig = plt.figure(figsize=(4.5, 4.5))
    ax = fig.add_subplot(111)
    plt.scatter(xis, yis, color=COLOR1, s=0.1)
    ax.set_aspect('equal', adjustable='box')
    ax.yaxis.set_tick_params(direction='in', which='both')
    ax.xaxis.set_tick_params(direction='in', which='both')

    plt.xlabel("X [kpc]")
    plt.ylabel("Y [kpc]")
    plt.title("Initial Particle Distribution")

    return print("go image/ check out " + plotName), plt.savefig('images/'+ plotName, dpi=200)

def plotPdfCdf():
    """
    Plot and save the image which contains the probability density function
    and the cumulative density function curves in directory 'images/'

    """
    plotName = 'E_iniColDenPdfCdf.png'
    random_samples = np.array([random.random() for i in range(5000)])
    samples = [invCDF(rand, SCALE, RM) for rand in random_samples]

    (figure, axes) = plt.subplots()
    xs = np.arange(0., RM, 0.01)
    yes = den(xs, SCALE)
    yps = pDF(xs, SCALE, RM)
    yns = cDF(xs, SCALE, RM)

    ys = [invHelper(i, 0.4, 3, 15) for i in xs]

    plt.plot(xs, yes, label="Column density", color=COLOR1, linewidth=2.5)
    plt.plot(xs, yps, label="PDF", color=COLOR2, linewidth=2.5)
    plt.plot(xs, yns, label="CDF ", color=COLOR3, linewidth=2.5)
    plt.plot(xs, ys, label="CDF monotonicity ", linestyle='--', color=COLOR4, linewidth=2.5, )
    plt.ylim(10 ** (-4), 5)
    plt.yscale("log")
    axes.yaxis.set_tick_params(direction='in', which='both')
    axes.xaxis.set_tick_params(direction='in', which='both')

    num_bins = 30
    n, bins, patches = plt.hist(samples, num_bins, facecolor=COLOR5, density=True, alpha=0.5)

    plt.title("Exponential Initial Column Surface Density \n scale = " + str(SCALE) + "Gyr")
    plt.xlabel("R [kpc]")
    plt.ylabel("Value")
    plt.legend()


    return print("go images/ to check " + plotName), plt.savefig('images/'+ plotName, dpi=200)
