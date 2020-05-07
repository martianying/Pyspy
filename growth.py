#--------------------GROWTH FUNCTION------------------#
import numpy as np
from phycon import *
import math as m

#----halo 
def haloGrow(ti, duration):
    ti = ti * gyrTos
    duration = duration * gyrTos
    
    if ti <= duration:
        td = 2.* ti / duration - 1.
        mag = 1. * (3. / 16. * td**5 - 5. / 8. * td**3 + 15. / 16. * td + 1./2.)
    else:
        mag = 1.
    return mag
    
    
#----cgArm 
def armEvol(ti, tpeakf = TPEAK / gyrTos, sigma = 1):
    tpeakf = tpeakf * gyrTos
    sigma = sigma * gyrTos
    ti = ti * gyrTos
    mag = np.exp( - (ti - tpeakf)**2 / (2. * sigma**2) )
    return mag
    
    
#----bar
def barGrowth(ti):
    taub = 6 * m.pi / PHIBAR
    t0 = 4 * taub   #time used to reach the half maximum of the total bar strength
    ti = ti * gyrTos
    mag = 1/2 + 1/2 * m.tanh((ti- t0)/taub)
    return mag
    

#----tidal
def tidalTerm(ti):
    ti = ti * gyrTos
    taut = 1 / PHIBAR
    t0 = 2.5 * taut
    mag = m.exp(-(ti - 5 * t0)**2 / (2 * (taut)**2))
    return mag
