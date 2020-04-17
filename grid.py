
# ---------------------- LOAD PACKAGES ---------------------#

import numpy as np
import math as m
import matplotlib.pyplot as plt
from matplotlib.pyplot import figure
from numpy import loadtxt
import seaborn as sns; sns.set()
import os

# ------------------ SET WORKING DIRECTORY -----------------#

os.chdir('/Users/veronicaplanck/Desktop/dataset/original')

# ------------------ SET UP THE PARAMETERS -----------------#

Rm = 18                                                  # Radius of the disc

reso_R = 80
dis_R = np.linspace(0, Rm, reso_R + 1)                       # Divide R into pieces

reso_theta = 60                                          # Divide θ into pieces
dis_theta = np.linspace(0, m.pi /2, reso_theta + 1)

raw_data_file = "blank80.txt"     # Path of raw data file which stores the x-y position                                                                information of each particle
# mass_data_file = "../mass/mass_" + raw_data_filabse
# den_data_file = "../density/den_" + raw_data_filabse
mass_data_file = "../mass/80mass.txt"                                # File that stores the particle numbers in each grid
den_data_file = "../density/80den.txt"                                # File that stores the mass_r density in each grid

# ---------------------READ INPUT DATA---------------------#
# Note:
# Extract x-y of each particle in raw_data_f and put into npos.
# npos = np.array([[x1,y1],[x2,y2],[x3,y3],...])

f = open(raw_data_file ,'r')
x = loadtxt(raw_data_file ,unpack=True ,usecols=[0])
y = loadtxt(raw_data_file ,unpack=True ,usecols=[1])
npos = np.array([[x[i] ,y[i]] for i in range(len(x))])


# ----------------------- FUNCTIONS -----------------------#

def ringDisc(array_of_pos, r_loc = None):
    """
    ringDisc(position_data: array_of_pos, int: r_loc)
    r_loc has to be within the range of the reso_R.

    This function READS: the position data and the id of target ring.

    When r_loc is not given:
    This function devides the whole disc into (reso_R -1) rings and RETURNS particle numbers
    in each ring.

    When r_loc => range(0, reso_r, 1):
    This function devides the whole disc into (reso_R -1) rings and RETURNS particle numbers
    in the selected ring and the detailed positions of those particles.

    """
    data_new = np.copy(array_of_pos)
    k = 0
    # Create an empty list to gather indexs of particles.
    data_basket = []

    while k < reso_R:
        # Find indexs of particles located between R[k] and R[k+1].
        np_in_a_ring = np.where( np.logical_and( np.sqrt( data_new[: ,0 ]**2 + data_new[: ,1 ]**2) <= dis_R[ k +1],
                                                 np.sqrt( data_new[: ,0 ]**2 + data_new[: ,1 ]**2 ) > dis_R[k] ) )

        # Put the above postions of particles as a list into the basket.
        data_basket = data_basket + [data_new[np_in_a_ring].tolist()]

        # Delete particles located between R[k] and R[k+1] to save calculation.
        data_new = np.delete(data_new, np_in_a_ring, 0)

        # Move on to the next outer ring.
        k = k + 1

    # return particle numbers in each rings.
    if r_loc is None :
        result = [len(data_basket[j]) for j in range(len(data_basket))]

    # return particle numbers in the selected ring and the detailed positions of those particles.
    else :
        result = len((data_basket[r_loc])), data_basket[r_loc]

    return result


def toQuadrant(array_of_pos, quadrant):
    """
    toQuadrant(position_data: array_of_pos, int: quadrant)
    quadrant has to be one of 1, 2, 3 and 4.

    This function READS: the position data and desired quadrant.

    This function divides the disc into the four quadrants and RETURNS particle number
    in the desired quadrant and the detailed particle postions data array_of_pos in that quadrant.

    """

    data_new = np.copy(array_of_pos)

    # Find the particles in the 1st quadrant.
    if quadrant == 0:
        particles_in_this_quadrant = np.where(np.logical_and(data_new[: ,1] > 0,
                                                             data_new[: ,0] > 0))

    # Find the particles in the 2nd quadrant.
    elif quadrant == 1:
        particles_in_this_quadrant = np.where(np.logical_and(data_new[: ,1] > 0,
                                                             data_new[: ,0] < 0))

    # Find the particles in the 3rd quadrant.
    elif quadrant == 2:
        particles_in_this_quadrant = np.where(np.logical_and(data_new[: ,1] < 0,
                                                             data_new[: ,0] < 0))

    # Find the particles in the 4th quadrant.
    else:
        particles_in_this_quadrant = np.where(np.logical_and(data_new[: ,1] < 0,
                                                             data_new[: ,0] > 0))

    return len(particles_in_this_quadrant[0]), array_of_pos[particles_in_this_quadrant]





def angleDisc(array_of_pos, theta_loc = None):
    """
    angleDisc(position_data: array_of_pos, int: r_loc)
    r_loc has to be within the range of the reso_theta.

    This function READS: the position data and the id of target ring.

    When theta_loc is not given:
    This function devides the whole disc into (reso_theta -1) pieces outwards and RETURNS particle numbers
    in each piece.

    When theta_loc => range(0, reso_theta, 1):
    This function devides the whole disc into (reso_theta -1) pieces outwards and RETURNS particle numbers
    in each piece and the detailed positions of those particles.

    """
    data_new = np.copy(array_of_pos)
    k = 0
    # Create an empty list to gather indexs of particles.
    data_basket = []  # here it's a tuple random

    while k < reso_theta:
        # Find indexs of particles located between theta[k] and theta[k+1].
        np_in_an_angle = np.where(np.logical_and(np.abs(np.arctan(data_new[: ,1 ] /data_new[: ,0])) <= dis_theta[ k +1],
                                                 np.abs(np.arctan(data_new[: ,1 ] /data_new[: ,0])) > dis_theta[k]))

        # Put the above postions of particles as a list into the basket.
        data_basket = data_basket + [data_new[np_in_an_angle].tolist()]

        # Delete particles located between theta[k] and theta[k+1] to save calculation.
        data_new = np.delete(data_new, np_in_an_angle, 0)

        # Move on to the next outer piece.
        k = k + 1

    # return particle numbers in each piece.
    if theta_loc is None:
        result = [len(data_basket[j]) for j in range(len(data_basket))]

    # return particle numbers in the selected piece and the detailed positions of those particles.
    else:
        result = len((data_basket[theta_loc])), data_basket[theta_loc]

    return result


def secterArea(loc):
    """
    secterArea(int: loc)
    loc shoud be in range(0, reso_theta, 1)

    This function READS the location of sector and RETURNS the area.

    """
    # Height of each sector
    height = reso_R
    # Outer base of the selected sector
    length_up = dis_theta[1] * dis_R[loc + 1]
    # Inter base of the selected sector
    length_down = dis_theta[1] * dis_R[loc]

    # Area of the selected sector
    result = 1/ 2 * (length_up + length_down) * height

    return result


def gridR(array_of_pos, quadrant, r_loc, theta_loc):
    """
    gridR(array_of_pos, quadrant, r_loc, theta_loc)
    theta_loc shoud be in range(0, reso_theta, 1)
    r_loc shoud be in range(0, reso_r, 1)

    This function READS: positions of particles, quadrant and the radial id and spatial id.
    It grids the position data into cells in the given quadrant.

    RETURNS: the number of particles in one cell (theta, R).

    """

    # Get all the particels in given quadrant.
    quatered_data = toQuadrant(array_of_pos, quadrant)[1]
    # Get the position data of particles in target radial id.
    ringed_data = ringDisc(quatered_data, r_loc)[1]
    # Get the number of particles in target radial id and spatial id.
    particle_number_in_cell = angleDisc(ringed_data, theta_loc)[0]

    return particle_number_in_cell


# start to write the density file
den_w = open(den_data_file, "w")

for i in range(4):
    for j in range(len(dis_R) - 1):
        # Ring particles in quadrant i at within angle j and get the postions.
        mass_within_R = angleDisc(ringDisc(toQuadrant(npos, i)[1], j)[1])
        den_of_cell_within_R = np.array(mass_within_R) / secterArea(j)
        # Write the mass_within_theta into file mass_data_file.
        for k in range(len(mass_within_R)):
            den_w.write("%f\t" % (den_of_cell_within_R[k]))

        # Start a new line to read particles in quadrant i at within angle j+1.
        den_w.write("\n")
den_w.close()


f = open("../density/80den.txt",'r')
lines = f.readlines()
loaded_data = []
for line in lines:
    loaded_data += [line.split()]
# close the file after reading the lines.
f.close()

def denThetaRadius(inarray):
    inarray = [inarray[ : reso_R], inarray[reso_R : 2*reso_R:], inarray[2 * reso_R : 3 * reso_R:], inarray[3 * reso_R :]]
    out_arr = inarray[0]
    for i in range(1, 4):
        out_arr = np.column_stack((out_arr, inarray[i]))
    return out_arr


[line.reverse() for line in loaded_data[reso_R : 2 * reso_R]];
[line.reverse() for line in loaded_data[3 * reso_R: ]];
read_den_data_reversed = np.array(loaded_data).astype(float)

sns.heatmap(denThetaRadius(read_den_data_reversed), cmap="Blues")









