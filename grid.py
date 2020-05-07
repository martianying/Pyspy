
# ---------------------- LOAD PACKAGES ---------------------#

import numpy as np
import math as m
import matplotlib.pyplot as plt
from matplotlib.pyplot import figure
from numpy import loadtxt
import ntpath
from phycon import *
import os

# ------------------ SET WORKING DIRECTORY -----------------#

os.chdir('/Users/veronicachang/Desktop/pyThesis')


# ----------------------- FUNCTIONS -----------------------#

def ringDisc(data, r_loc):
    data_new = np.copy(data)
    k = 0
    data_basket = []

    while k < RESO_R:
        np_in_a_ring = np.where(np.logical_and(np.sqrt(data_new[: ,0 ]**2 + data_new[: ,1 ]**2) <= DIS_R[ k +1],
                                               np.sqrt(data_new[: ,0 ]**2 + data_new[: ,1 ]**2 ) >DIS_R[k]))

        data_basket = data_basket + [data_new[np_in_a_ring].tolist()]
        data_new = np.delete(data_new, np_in_a_ring, 0)

        k = k + 1
    if r_loc >= 0:
        result = len((data_basket[r_loc])), data_basket[r_loc]
    else:
        result = [len(data_basket[j]) for j in range(len(data_basket))]
    return result

#   len(data_basket[n_radii])

def toQuadrant(data, nq):
    data_new = np.copy(data)
    if nq == 0:
        test = np.where(np.logical_and(data_new[: ,1] > 0,
                                       data_new[: ,0] > 0))
    elif nq == 1:
        test = np.where(np.logical_and(data_new[: ,1] > 0,
                                       data_new[: ,0] < 0))
    elif nq == 2:
        test = np.where(np.logical_and(data_new[: ,1] < 0,
                                       data_new[: ,0] < 0))
    else:
        test = np.where(np.logical_and(data_new[: ,1] < 0,
                                       data_new[: ,0] > 0))
    return len(test[0]) ,data[test]



def sliceDisc(data, theta_loc):
    data_new = np.copy(data)
    k = 0
    data_basket = []  # here it's a tuple random

    while k < RESO_THETA:
        np_in_an_angle = np.where(np.logical_and(np.abs(np.arctan(data_new[: ,1 ] /data_new[: ,0])) <= DIS_THETA[ k +1],
                                                 np.abs(np.arctan(data_new[: ,1 ] /data_new[: ,0])) > DIS_THETA[k]))
        data_basket = data_basket + [data_new[np_in_an_angle].tolist()]
        k = k + 1
    if theta_loc >= 0:
        result = len((data_basket[theta_loc])), data_basket[theta_loc]
    else:
        result = [len(data_basket[j]) for j in range(len(data_basket))]
    return result


def secterArea(loc):
    width = RM / RESO_R
    length_up = DIS_THETA[1] * DIS_R[loc +1]
    length_down = DIS_THETA[1] * DIS_R[loc]
    result = 1/ 2 * (length_up + length_down) * width
    return result


def gridR(data, nq, r_loc, theta_loc):
    quatered_data = toQuadrant(data, nq)[1]
    ringed_data = ringDisc(quatered_data, r_loc)[1]
    cell_num = sliceDisc(ringed_data, theta_loc)[0]
    return cell_num  # return the number of particles in one cell (theta, R).


def getDenFileOf(raw_data_file):
    raw_data_file_base = ntpath.basename(raw_data_file)
    x = loadtxt(raw_data_file, comments='#', unpack=True, usecols=[0])
    y = loadtxt(raw_data_file, comments='#', unpack=True, usecols=[1])
    npos = np.array([[x[i], y[i]] for i in range(len(x))])

    np_data_file = "mass/mass_" + raw_data_file_base  # File that stores the particle numbers in each grid
    den_data_file = "den/den_" + raw_data_file_base  # File that stores the mass_r density in each grid
    mass_w = open(np_data_file, "w")

    for i in range(4):
        for j in range(len(DIS_THETA) - 1):
            np_within_theta = ringDisc(sliceDisc(toQuadrant(npos, i)[1], j)[1], -1)
            for k in range(len(np_within_theta)):
                mass_w.write("%i\t" % (np_within_theta[k]))
            mass_w.write("\n")
    mass_w.close()

    mass_r = open(np_data_file, "r")
    den_w = open(den_data_file, "w")

    for i in range(len(DIS_R) - 1):
        test = loadtxt(np_data_file, unpack=True, usecols=[i])
        den_for_this_col = test / secterArea(i)
        for k in range(len(den_for_this_col)):
            den_w.write("%f\t" % (den_for_this_col[k]))
        den_w.write("\n")
    mass_r.close()
    den_w.close()

    print("Procession for: " + raw_data_file + " is completed.")


def adjustFile(denfile_name, CONTRAST='y'):
    g = open(denfile_name, 'r')
    lines = g.readlines()
    loaded_data = []
    for line in lines:
        sline = line.split()

        rev2 = sline[RESO_THETA:2 * RESO_THETA]
        rev2.reverse()
        sline[RESO_THETA:2 * RESO_THETA] = rev2

        rev3 = sline[3 * RESO_THETA:]
        rev3.reverse()
        sline[3 * RESO_THETA:] = rev3

        loaded_data += [sline]

    # close the file after reading the lines.
    g.close()
    if CONTRAST == 'n' or CONTRAST == 'N':
        dendata = np.array(loaded_data).astype(float)
    else:
        dendata = np.array(loaded_data).astype(float)
        for i in range(len(dendata)):
            dendata[i] /= np.mean(dendata)
    return dendata


def plotHeatmap(denFileName):
    import seaborn as sns
    #here the adjustFile function has already directed the directory.
    adjusted_denFile = adjustFile(denFileName,'n')
    sns.heatmap(adjusted_denFile, cmap="Blues")
    plt.title("Heatmap of grided disc(To be customized)")
    plt.xlabel("DIS_THETA")
    plt.ylabel("DIS_R")
    return plt.show()


# plot the surface density along the radius at theta = 0.
def denAlongR(denfile_path_name, scolor, lab):
    # figure(dpi=300, facecolor='w', edgecolor='k')
    f = open(denfile_path_name)
    denRi = loadtxt(f, unpack=True, usecols=[0])
    radius = DIS_R[:-1]

    # plt.plot(radius, denRi * 0.01, linewidth=2.5)
    plt.scatter(radius, denRi * 0.01, color=scolor, s=25, label=lab)
    plt.yscale("log")
    #sns.set(style="ticks", palette="muted")
    plt.xlabel('R [kpc]')
    plt.ylabel(r'Column Density  [g/$kpc^{2}$]')
    # plt.ylim(0, y_lim)

    plt.title('Initial Radial Column Density')
    plt.legend()

    return plt.savefig(denfile_path_name[:-4] + '.png', dpi=200)




