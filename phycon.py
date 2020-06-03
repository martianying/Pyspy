#-----------------PHYSICAL CONSTANTS-----------------#
import numpy as np
import math as m

#------To grid the disc
RM = 18                                                  # Radius of the disc 

RESO_R = 80      
DIS_R = np.linspace(0, RM, RESO_R + 1)                       # Divide R into pieces

RESO_THETA = 60                                          # Divide θ into pieces        
DIS_THETA = np.linspace(0, m.pi/2, RESO_THETA + 1)



#------Units Conversion
kpcTom = 3.086*10**19           #distance 
cmTom = 10**(-2)
kmTom = 10**3

gTokg = 10**3                   #weight 
solarM = 1.98*10**30            #solar mass
Gcode =  6.7 * 10**(-11)        #gravitational constant

PARNUM = 10**6
PARTICLE_MASS = 10**10 / PARNUM   #mass of gas particle in kg

gyrTos = 3.15 * 10**16          #time 
T_UNIT_PHANTOM = 0.0468         #Gyr

#------Halo Constants
MHALO = 6.4 * 10**10 * solarM
OMEH = 20 * kmTom / kpcTom

RHALOMAX = 12. * kpcTom
RHALO = 6. * kpcTom
HALOCI  = 1. / ( RHALOMAX / RHALO - m.atan( RHALOMAX / RHALO ))

#------Disc Constants
MDISC = 8.56 * 10**10 * solarM
ADISC = 5.32 * kpcTom
BDISC = 0.25 * kpcTom

#------Bulge Constants
MBULGE = 1.4 * 10**10 * solarM
BBULGE = 0.39 * kpcTom

#------Bar Constants
ALPHA2 = 0.7190                  #constants for LMABAR = 2 * LMbbar = 2 * LMcbar
BETA2 = 0.6901
PHIBAR = 11 * kmTom / kpcTom     #rotational speed of bar
LMABAR = 4 * kpcTom              #characteristic length of bar
    
MBAR = 3. * 10**8 * solarM       #(*in kg*)

#------Tidal Constants
PHIBAR_T = 11 * kmTom / kpcTom
MBAR_T = 3.30 * 10**9 * solarM 

#------Arm Constants
RAMR0 = 8 * kpcTom
RARMS = 7 * kpcTom
RHOARM = 2.1289 * 10**(-24) / (gTokg * cmTom**3)
TPEAK = 3 * gyrTos
NARM = 2
PITCHARM = m.pi / 180 * 15
HARM = 0.8 * kpcTom
CnARM = [8/(3 * m.pi), 1/2, 8/(15 * m.pi)]
CGARM = -4 * m.pi * Gcode * HARM * RHOARM
