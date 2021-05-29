#----------------------------------------------------------------------------#
# MODULE: roCurve                                                            #
#                                                                            #
# CONTENT:                                                                   #
#        gravitational potential components of the gaseous disc              #
#     and their contribution on the velocities of gas particles.             #
#                                                                            #
# AIM:                                                                       #
#      1) plot the predicted rotation curve of the disc with/without bulge   #
#      2) read PHANTOM splash files and plot the actual rotation curve in    #
#      your simulation, so as to do the model setup check.                   #
#                                                                            #
# DEPENDENCIES: phycon                                                       #
#                                                                            #
#----------------------------------------------------------------------------#


import matplotlib.pyplot as plt
from src.phycon import *
from numpy import loadtxt
import os

# ~~~~~~~~~~~~~~~~~~DISC~~~~~~~~~~~~~~~~~~~~~~
def disc(r_kpc):
    """
    Calculate the velocity at a given radius in the disc due to
    the disc potential.

    Argument:
    r_kpc: radius in kpc

    Returns:
    squared velocity at r_kpc in km
    """
    
    r_m = r_kpc * kpcTom 
    Up = Gcode * MDISC * r_m**2
    Down = (ADISC + BDISC)**2 + r_m**2
    vel2_in_m = Up / Down**(3/2)
    vel2_in_km = vel2_in_m / kmTom**2
    
    return vel2_in_km

# ~~~~~~~~~~~~~~~~BULGE~~~~~~~~~~~~~~~~~~~~~
def bulge(r_kpc):
    """
    Calculate the velocity at a given radius in the disc due to
    the bulge potential.

    Argument:
    r_kpc: radius in kpc

    Returns:
    squared velocity at r_kpc in km
    """
    
    r_m = r_kpc * kpcTom 
    Up = Gcode * MBULGE * r_m**2
    Down = r_m**2 + BBULGE**2
    vel2_in_m = Up / Down**(3 / 2)
    vel2_in_km = vel2_in_m / kmTom**2
    
    return vel2_in_km
    

# ~~~~~~~~~~~~~~~~~HALO~~~~~~~~~~~~~~~~~~~~~~
def halo(r_kpc):
    """
    Calculate the velocity at a given radius in the disc due to the
    halo potential

    Argument:
    r_kpc: radius in kpc

    Returns:
    squared velocity at r_kpc in km
    """

    r_m = r_kpc * 3.086 * 10**19
    Up = MHALO * Gcode * ( np.arctan(r_m / RHALO) - r_m / RHALO)
    Down = r_m * (RHALOMAX / RHALO - np.arctan(RHALOMAX / RHALO))
    vel2_in_m = - Up / Down
    vel2_in_km = vel2_in_m / kmTom**2
    
    
    return vel2_in_km
    

# ~~~~~~~~~~~~~~ROTATION VELOCITY~~~~~~~~~~~~~~
def rotation(r_kpc):
    """
    Calculate the velocity of a gas particle at a given radius in the disc.

    Argument:
    r_kpc: radius of the disc in kpc

    Returns:
    squared rotational velocity with bulge-free at r_kpc in km
    """

    Vc_2_in_km = disc(r_kpc) + halo(r_kpc)

    return Vc_2_in_km

def addbulge(r_kpc):
    """
    Calculate the velocity at a given radius in the disc with bulge
    taken into account.

    Argument:
    r_kpc: radius of the disc in kpc

    Returns:
    squared rotational velocity with bulge at r_kpc in km
    """
    Vc_2_in_km = disc(r_kpc) + halo(r_kpc) + bulge(r_kpc)

    return Vc_2_in_km


#~~~~~~~~~~~~~~~~PLOT~~~~~~~~~~~~~~~~~~~~~~~~
def rotationCurve(funcName, COLOR, LABEL,  STYLE, filename = None, ALPHA = 0.65):
    """
    Calculate the velocity at a given radius in the disc.

    Argument:
    funcName: type of calculated squared velocity (eg. disc, bulge, halo, rotation, etc.)
    COLOR:    color scheme (better use the one you defined in phycon)
    LABEL:    label of the curve
    STYLE:    style of the curve
    filename: txt file containing gas particles info obtained from your simulations results
    calculated by PHANTOM, converted by the SPLASH
    ALPHA:    transparency of the curve, for aesthetic reason

    Returns:
    rotation curve plot. (radius[kpc] vs velocity[km/s])
    """

    # scatter plot of simulation results
    if filename != None:
        unit_p = 1
        unit_v = 10**5
        parnum = 2000
        os.chdir('/Users/veronicachang/Desktop/pyThesis/datasets')

        raw_data_file = filename

        # load the positions and velocities of gas particles in the given filename
        xi, yi, zi = [loadtxt(raw_data_file, comments = '#', unpack = True, usecols = [i]) for i in range(3)]
        vxi, vyi, vzi = [loadtxt(raw_data_file, comments = '#', unpack= True, usecols= [i]) for i in range(6, 9)]

        # only show parnum particles in the rotation curve for time saving
        xiSub, yiSub, ziSub = xi[:parnum], yi[:parnum], zi[:parnum]
        vxiSub, vyiSub, vziSub = vxi[:parnum], vyi[:parnum], vzi[:parnum]

        # unit conversion :  from PHANTOM unit to km/s and kpc
        rlist = list([xiSub, yiSub, ziSub])
        vrlist = list([vxiSub, vyiSub, vziSub])
        rlist_kpc = [np.array(ele) * unit_p for ele in rlist]
        vrlist_kms = [np.array(ele) / unit_v for ele in vrlist]

        # calculte the distance to the origin: Rs and the velocities: Vs, of gas particles
        Rs = np.sqrt(rlist_kpc[0] ** 2 + rlist_kpc[1] ** 2 + rlist_kpc[2] ** 2)
        Vs = np.sqrt(vrlist_kms[0] ** 2 + vrlist_kms[1] ** 2 + vrlist_kms[2] ** 2)

        # scatter plot each particle on the radius[kpc] vs velocity[km/s] curve
        rotationPlot = plt.scatter(Rs, Vs, label = 'Simulation Samples', s = 40, color = "#A569BD")

    # plot the predicted rotation curve
    Vci = []

    r_kpc = np.arange(0.001, RM, 0.01)

    for i in range(len(r_kpc)):
        Vci = Vci + [funcName(r_kpc[i])]
    Vci = np.sqrt(np.array(Vci))

    plt.plot(r_kpc, Vci, color=COLOR, linewidth=2.8, label=LABEL, linestyle = STYLE, alpha = ALPHA)
    plt.xlabel("Radius [kpc]")
    plt.ylabel("Velocity [km/s]")
    plt.legend(loc='lower right')
    plt.title("Rotation Curve")

    return plt.savefig('../images/'+ 'Rc.png', dpi=200)
