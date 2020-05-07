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
    Up = Gcode * MDISC * r_m**2
    Down = (ADISC + BDISC)**2 + r_m**2
    vel2_in_m = Up / Down**(3/2)
    vel2_in_km = vel2_in_m / kmTom**2
    
    return vel2_in_km

#------BULGE
def bulge(r_kpc):
    
    r_m = r_kpc * kpcTom 
    Up = Gcode * MBULGE * r_m**2
    Down = r_m**2 + BBULGE**2
    vel2_in_m = Up / Down**(3 / 2)
    vel2_in_km = vel2_in_m / kmTom**2
    
    return vel2_in_km
    

#------HALO
def halo(r_kpc):
    
    r_m = r_kpc * 3.086 * 10**19
    Up = MHALO * Gcode * ( np.arctan(r_m / RHALO) - r_m / RHALO)
    Down = r_m * (RHALOMAX / RHALO - np.arctan(RHALOMAX / RHALO))
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
def rotationCurve(funcName, COLOR, LABEL,  STYLE, filename = None, ALPHA = 0.65):

    if filename != None:
        unit_p = 1
        unit_v = 10**5
        parnum = 2000
        os.chdir('/Users/veronicachang/Desktop/pyThesis/datasets')

        raw_data_file = filename

        xi, yi, zi = [loadtxt(raw_data_file, comments = '#', unpack = True, usecols = [i]) for i in range(3)]
        vxi, vyi, vzi = [loadtxt(raw_data_file, comments = '#', unpack= True, usecols= [i]) for i in range(6, 9)]

        xiSub, yiSub, ziSub = xi[:parnum], yi[:parnum], zi[:parnum]
        vxiSub, vyiSub, vziSub = vxi[:parnum], vyi[:parnum], vzi[:parnum]

        rlist = list([xiSub, yiSub, ziSub])
        vrlist = list([vxiSub, vyiSub, vziSub])
        rlist_kpc = [np.array(ele) * unit_p for ele in rlist]
        vrlist_kms = [np.array(ele) / unit_v for ele in vrlist]

        Rs = np.sqrt(rlist_kpc[0] ** 2 + rlist_kpc[1] ** 2 + rlist_kpc[2] ** 2)
        Vs = np.sqrt(vrlist_kms[0] ** 2 + vrlist_kms[1] ** 2 + vrlist_kms[2] ** 2)

        rotationPlot = plt.scatter(Rs, Vs, label = 'Simulation Samples', s = 40, color = "#A569BD")

    Vci = []

    r_lim = 20
    r_kpc = np.arange(0.001, r_lim, 0.01)

    for i in range(len(r_kpc)):
        Vci = Vci + [funcName(r_kpc[i])]
    Vci = np.sqrt(np.array(Vci))

    plt.plot(r_kpc, Vci, color=COLOR, linewidth=2.8, label=LABEL, linestyle = STYLE, alpha = ALPHA)
    plt.xlabel("Radius [kpc]")
    plt.ylabel("Velocity [km/s]")
    plt.legend(loc='lower right')
    plt.title("Rotation Curve")

    return plt.savefig('../images/'+ 'Rc.png', dpi=200)







