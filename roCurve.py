#------------------ROTATION CURVE--------------------#

import numpy as np
import matplotlib.pyplot as plt
from phycon import *
from numpy import loadtxt
import os

#------DISC
def disc(r_kpc):
    '''
    IN:radius in kpc; OUT: velocity in km/s
    '''
    
    r_m = r_kpc * kpcTom 
    Up = Gcode * Md * r_m**2
    Down = (ad + bd)**2 + r_m**2
    vel2_in_m = Up / Down**(3/2)
    vel2_in_km = vel2_in_m / kmTom**2
    
    return vel2_in_km

#------BULGE
def bulge(r_kpc):
    
    r_m = r_kpc * kpcTom 
    Up = Gcode * Mb * r_m**2
    Down = r_m**2 + bb**2
    vel2_in_m = Up / Down**(3 / 2)
    vel2_in_km = vel2_in_m / kmTom**2
    
    return vel2_in_km
    

#------HALO
def halo(r_kpc):
    
    r_m = r_kpc * 3.086 * 10**19
    Up = Mhalo * Gcode * ( np.arctan(r_m / rhalo) - r_m / rhalo)
    Down = r_m * (rhaloMax / rhalo - np.arctan(rhaloMax / rhalo))
    vel2_in_m = - Up / Down
    vel2_in_km = vel2_in_m / kmTom**2
    
    
    return vel2_in_km
    

#------TOTAL---#
def rotation(r_kpc):
    Vc_2_in_km = disc(r_kpc) + halo(r_kpc)
    return Vc_2_in_km

def addbulge(r_kpc):
    Vc_2_in_km = disc(r_kpc) + halo(r_kpc) + bulge(r_kpc)
    return Vc_2_in_km


#-------the plot 
#rotationCurve(rotation, black, "Rotation Curve", '+', 'blank80')
def rotationCurve(function, COLOR, LABEL, STYLE, filename = None):

    if filename != None:
        unit_p = 1
        unit_v = 10**5
        parnum = 2000
        os.chdir('/Users/veronicaplanck/Desktop/pyAnalysis/datafiles')

        raw_data_file = str(filename) + '.txt'

        xi, yi, zi = [loadtxt(raw_data_file, unpack = True, usecols = [i]) for i in range(3)]
        vxi, vyi, vzi = [loadtxt(raw_data_file, unpack= True, usecols= [i]) for i in range(6, 9)]

        xiSub, yiSub, ziSub = xi[:parnum], yi[:parnum], zi[:parnum]
        vxiSub, vyiSub, vziSub = vxi[:parnum], vyi[:parnum], vzi[:parnum]

        rlist = list([xiSub, yiSub, ziSub])
        vrlist = list([vxiSub, vyiSub, vziSub])
        rlist_kpc = [np.array(ele) * unit_p for ele in rlist]
        vrlist_kms = [np.array(ele) / unit_v for ele in vrlist]

        Rs = np.sqrt(rlist_kpc[0] ** 2 + rlist_kpc[1] ** 2 + rlist_kpc[2] ** 2)
        Vs = np.sqrt(vrlist_kms[0] ** 2 + vrlist_kms[1] ** 2 + vrlist_kms[2] ** 2)

        rotationPlot = plt.scatter(Rs, Vs)

    Vci = []

    r_lim = 20
    r_kpc = np.arange(0.001, r_lim, 0.01)

    for i in range(len(r_kpc)):
        Vci = Vci + [function(r_kpc[i])]
    Vci = np.sqrt(np.array(Vci))

    plt.plot(r_kpc, Vci, STYLE, alpha=0.8, color=COLOR, linewidth=3, label=LABEL)
    plt.xlabel("radius_in_kpc")
    plt.ylabel("velocity_in_km/s")
    plt.legend(loc='lower right')
    plt.title("ROTATION CURVE")

    return plt.show()








