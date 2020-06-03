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
    #return the quantity in m^2.
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


def plotHeatmap(denFileName, TITLE):
    import seaborn as sns

    # here the adjustFile function has already directed the directory.
    adjusted_denFile = adjustFile(denFileName, 'n')
    # fig, ax = plt.subplots(figsize=(7,5))
    heat_map = sns.heatmap(adjusted_denFile, cmap="RdPu",
                           cbar_kws={'label': 'Dnesity Contrast'}, vmin=0, vmax=1700)

    cbar = heat_map.collections[0].colorbar
    cbar.set_ticks([0, 1700 * 1 / 4, 1700 * 2 / 4, 1700 * 3 / 4, 1700])
    cbar.set_ticklabels(['0', '25%', '50%', '75%', '100%'])
    cbar.ax.tick_params(labelsize=8, grid_alpha=0.5, direction='in')

    plt.xticks([0, RESO_THETA, 2 * RESO_THETA, 3 * RESO_THETA, 4 * RESO_THETA],
               ['0', '$\pi$ / 2', '$\pi$', '$3\pi / 2$', "2$\pi$"], rotation=0)

    plt.yticks([0, RESO_R / RM * 5, RESO_R / RM * 10, RESO_R / RM * 15],
               ['0', '5', '10', '15'], rotation=0)

    plt.tick_params(direction="in")
    plt.title(TITLE, x=0.5, y=1.05)
    plt.xlabel(r"$\theta$ [rad]")
    plt.ylabel("R [kpc]")

    return plt.savefig('images/' + TITLE, dpi=200), print("finished")



# plot the surface density along the radius at theta = 0.
def denAlongR(denfile_path_name_list, scolor_list, lab_list):
    (figure, axes) = plt.subplots()
    for i in range(len(denfile_path_name_list)):
        f = open(denfile_path_name_list[i])
        denR_num  = loadtxt(f, unpack=True, usecols=[0])
        radius = DIS_R[:-1]
        # NumberDensity * ParticleMass (Mo/kpc^2)
        denRi = denR_num * PARTICLE_MASS

        plt.scatter(radius, denRi, color=scolor_list[i], s=25, label=lab_list[i])
        plt.yscale("log")

        plt.xlabel('R [kpc]')
        plt.ylabel(r'Column Density  [$M_{o}$ / $kpc^{2}$]')
        plt.tick_params(direction="in")

        axes.yaxis.set_tick_params(direction='in', which='both')
        plt.ylim(10, 10**10)

        plt.title('Initial Radial Column Density')
        plt.legend()
    return plt.savefig(denfile_path_name_list[0][:-4] + '.png', dpi=200), print("go den/ to check your output.")
