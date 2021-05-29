#---------------------------------------------------------------------------#
# MODULE: phycon                                                            #
#                                                                           #
# CONTENT:                                                                  #
#      1) parameters of griding operations on disc                          #
#      2) units conversion (eg. distance, time, mass)                       #
#      3) universal constants (eg. solar mass, gravitational constant)      #
#      4) parameters of galactic component (eg, disc, halo, bulge)          #
#      5) parameters of each analytical potential                           #
#         (eg. bar, arms, tidal component)                                  #
#                                                                           #
# AIM:                                                                      #
#      setup the physics environment                                        #
# DEPENDENCIES: None                                                        #
#                                                                           #
#---------------------------------------------------------------------------#

import numpy as np
import math as m

#---------------------------------------------------------------------
# COLOR SCHEME
# (In order to have a consistent color scheme and avoid redundancy,
# we set up our customized colors here.)
#---------------------------------------------------------------------
COLOR1 = "#990000"
COLOR2 = "#A569BD"
COLOR3 = "#004d99"
COLOR4 = "#D35400"
COLOR5 = "#D7BDE2"

#---------------------------------------------------------------------
# PARAMETERS OF HOW TO GRID THE DISC
#---------------------------------------------------------------------

RM = 18                                                      # Radius of the disc in kpc

RESO_R = 80      
DIS_R = np.linspace(0, RM, RESO_R + 1)                       # divide R into RESO_R pieces

RESO_THETA = 60                                              # divide θ into RESO_THETA pieces
DIS_THETA = np.linspace(0, m.pi/2, RESO_THETA + 1)


#---------------------------------------------------------------------
# UNITS CONVERSION
#---------------------------------------------------------------------
kpcTom = 3.086*10**19                                        # distance conversion    :  kpc to m
cmTom = 10**(-2)                                             # distance conversion    :  cm to m
kmTom = 10**3                                                # distance conversion    :  km to m

gTokg = 10**3                                                # weight conversion      :  g to kg
solarM = 1.98*10**30                                         # solar mass constant    :  kg
Gcode =  6.7 * 10**(-11)                                     # gravitational constant :  N*m^2/(kg)^2

PARNUM = 10**6                                               # number of particles in simulations in PHANTOM
PARTICLE_MASS = 10**10 / PARNUM                              # mass of gas particle   :  kg

gyrTos = 3.15 * 10**16                                       # time conversion        :  Gyr to second
T_UNIT_PHANTOM = 0.0468                                      # time unit in PHANTOM   :  Gyr


#---------------------------------------------------------------------
# PARAMETERS OF GALACTIC DISC COMPONENTS
# REFERENCES:

#~~~~~~~~~~~~~HALO~~~~~~~~~~~~~~~~~
MHALO = 6.4 * 10**10 * solarM
OMEH = 20 * kmTom / kpcTom

RHALOMAX = 12. * kpcTom
RHALO = 6. * kpcTom
HALOCI  = 1. / ( RHALOMAX / RHALO - m.atan( RHALOMAX / RHALO ))

#~~~~~~~~~~~~~DISC~~~~~~~~~~~~~~~~~
MDISC = 8.56 * 10**10 * solarM
ADISC = 5.32 * kpcTom
BDISC = 0.25 * kpcTom

#~~~~~~~~~~~~~BULGE~~~~~~~~~~~~~~~~
MBULGE = 1.4 * 10**10 * solarM
BBULGE = 0.39 * kpcTom

#---------------------------------------------------------------------
# PARAMETERS OF GALACTIC ANALYTICAL POTENTIAL
# REFERENCES:

#~~~~~~~~~~~~~~~BAR~~~~~~~~~~~~~~~~
ALPHA2 = 0.7190                                              # constants for LMABAR = 2 * LMbbar = 2 * LMcbar
BETA2 = 0.6901
PHIBAR = 11 * kmTom / kpcTom                                 # rotational speed of bar
LMABAR = 4 * kpcTom                                          # characteristic length of bar
    
MBAR = 3. * 10**8 * solarM                                   # unit  :  kg

#~~~~~~~~~~TIDAL COMPONENT~~~~~~~~~
PHIBAR_T = 11 * kmTom / kpcTom
MBAR_T = 3.30 * 10**9 * solarM 

#~~~~~~~~~~~~~~ARMS~~~~~~~~~~~~~~~~
RAMR0 = 8 * kpcTom
RARMS = 7 * kpcTom
RHOARM = 2.1289 * 10**(-24) / (gTokg * cmTom**3)
TPEAK = 3 * gyrTos
NARM = 2
PITCHARM = m.pi / 180 * 15
HARM = 0.8 * kpcTom
CnARM = [8/(3 * m.pi), 1/2, 8/(15 * m.pi)]
CGARM = -4 * m.pi * Gcode * HARM * RHOARM
