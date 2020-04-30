
# ---------------------- LOAD PACKAGES ---------------------#

import numpy as np
import math as m
import matplotlib.pyplot as plt
from matplotlib.pyplot import figure
from numpy import loadtxt

import os

Rm = 18                                                  # Radius of the disc

reso_R = 80
dis_R = np.linspace(0, Rm, reso_R + 1)                       # Divide R into pieces

reso_theta = 60                                          # Divide θ into pieces
dis_theta = np.linspace(0, m.pi /2, reso_theta + 1)

# ------------------ SET WORKING DIRECTORY -----------------#

os.chdir('/Users/veronicachang/Desktop/pyThesis')


# ----------------------- FUNCTIONS -----------------------#

def ringDisc(data, r_loc):
    data_new = np.copy(data)
    k = 0
    data_basket = []

    while k < reso_R:
        np_in_a_ring = np.where(np.logical_and(np.sqrt(data_new[: ,0 ]**2 + data_new[: ,1 ]**2) <= dis_R[ k +1],
                                               np.sqrt(data_new[: ,0 ]**2 + data_new[: ,1 ]**2 ) >dis_R[k]))

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

    while k < reso_theta:
        np_in_an_angle = np.where(np.logical_and(np.abs(np.arctan(data_new[: ,1 ] /data_new[: ,0])) <= dis_theta[ k +1],
                                                 np.abs(np.arctan(data_new[: ,1 ] /data_new[: ,0])) > dis_theta[k]))
        data_basket = data_basket + [data_new[np_in_an_angle].tolist()]
        k = k + 1
    if theta_loc >= 0:
        result = len((data_basket[theta_loc])), data_basket[theta_loc]
    else:
        result = [len(data_basket[j]) for j in range(len(data_basket))]
    return result


def secterArea(loc):
    width = Rm / reso_R
    length_up = dis_theta[1] * dis_R[loc +1]
    length_down = dis_theta[1] * dis_R[loc]
    result = 1/ 2 * (length_up + length_down) * width
    return result


def gridR(data, nq, r_loc, theta_loc):
    quatered_data = toQuadrant(data, nq)[1]
    ringed_data = ringDisc(quatered_data, r_loc)[1]
    cell_num = sliceDisc(ringed_data, theta_loc)[0]
    return cell_num  # return the number of particles in one cell (theta, R).


def getDenFileOf(raw_data_file):
    raw_data_file = "datasets/" + raw_data_file
    f = open(raw_data_file, 'r')
    x = loadtxt(raw_data_file, comments='#', unpack=True, usecols=[0])
    y = loadtxt(raw_data_file, comments='#', unpack=True, usecols=[1])
    npos = np.array([[x[i], y[i]] for i in range(len(x))])

    np_data_file = "mass/mass_" + raw_data_file  # File that stores the particle numbers in each grid
    den_data_file = "den/den_" + raw_data_file  # File that stores the mass_r density in each grid
    mass_w = open(np_data_file, "w")

    for i in range(4):
        for j in range(len(dis_theta) - 1):
            np_within_theta = ringDisc(sliceDisc(toQuadrant(npos, i)[1], j)[1], -1)
            for k in range(len(np_within_theta)):
                mass_w.write("%i\t" % (np_within_theta[k]))
            mass_w.write("\n")
    mass_w.close()

    mass_r = open(np_data_file, "r")
    den_w = open(den_data_file, "w")

    for i in range(len(dis_R) - 1):
        test = loadtxt(np_data_file, unpack=True, usecols=[i])
        den_for_this_col = test / secterArea(i)
        for k in range(len(den_for_this_col)):
            den_w.write("%f\t" % (den_for_this_col[k]))
        den_w.write("\n")
    mass_r.close()
    den_w.close()

    print("Process for" + raw_data_file + " is completed.")


# plot the surface density along the radius at theta = 0.
def denAlongR(denfile_path_name, scolor, lab):
    # figure(dpi=300, facecolor='w', edgecolor='k')
    f = open(denfile_path_name)
    denRi = loadtxt(f, unpack=True, usecols=[0])
    radius = dis_R[:-1]

    # plt.plot(radius, denRi * 0.01, linewidth=2.5)
    plt.scatter(radius, denRi * 0.01, color=scolor, s=25, label=lab)
    plt.yscale("log")
    #sns.set(style="ticks", palette="muted")
    plt.xlabel('R kpc')
    plt.ylabel(r'Column Density  g/$kpc^{2}$')
    # plt.ylim(0, y_lim)

    plt.title('Initial Radial Column Density')
    plt.legend()

    return plt.savefig(denfile_path_name[:-4] + '.png', dpi=200)


def denHeatmap(denfile_name):
    import seaborn as sns
    g = open(denfile_name, 'r')
    lines = g.readlines()
    loaded_data = []
    for line in lines:
        sline = line.split()

        rev2 = sline[reso_theta:2 * reso_theta]
        rev2.reverse()
        sline[reso_theta:2 * reso_theta] = rev2

        rev3 = sline[3 * reso_theta:]
        rev3.reverse()
        sline[3 * reso_theta:] = rev3

        loaded_data += [sline]

    # close the file after reading the lines.
    g.close()
    dendata = np.array(loaded_data).astype(float)
    return sns.heatmap(dendata, cmap="Blues")

