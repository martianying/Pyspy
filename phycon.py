#-----------------PHYSICAL CONSTANTS-----------------#
import numpy as np
import math as m

#------To grid the disc
Rm = 18                                                  # Radius of the disc 

reso_R = 80      
dis_R = np.linspace(0, Rm, reso_R + 1)                       # Divide R into pieces

reso_theta = 60                                          # Divide θ into pieces        
dis_theta = np.linspace(0, m.pi/2, reso_theta + 1)

#------Units Conversion
kpcTom = 3.086*10**19           #distance 
cmTom = 10**(-2)
kmTom = 10**3

gTokg = 10**3                   #weight 
solarM = 1.98*10**30            #solar mass
Gcode =  6.7 * 10**(-11)        #gravitational constant

gyrTos = 3.15 * 10**16          #time 
t_unit_phantom = 0.0468

#------Halo Constants
Mhalo = 6.4 * 10**10 * solarM
omeH = 20 * kmTom / kpcTom

rhaloMax = 12. * kpcTom
rhalo = 6. * kpcTom
HaloCi  = 1. / ( rhaloMax / rhalo - m.atan( rhaloMax / rhalo ))

#------Disc Constants
Md = 8.56 * 10**10 * solarM
ad = 5.32 * kpcTom
bd = 0.25 * kpcTom

#------Bulge Constants
Mb = 1.4 * 10**10 * solarM
bb = 0.39 * kpcTom

#------Bar Constants
alpha2 = 0.7190                  #constants for LMabar = 2 * LMbbar = 2 * LMcbar
beta2 = 0.6901
pi = m.pi
phibar = 11 * kmTom / kpcTom     #rotational speed of bar
LMabar = 4 * kpcTom              #characteristic length of bar
    
Mbar = 3. * 10**8 * solarM       #(*in kg*)

#------Tidal Constants
phibar_t = 11 * kmTom / kpcTom
Mbar_t = 3.30 * 10**9 * solarM 
