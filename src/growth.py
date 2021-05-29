#----------------------------------------------------------------------------#
# MODULE: roCurve                                                            #
#                                                                            #
# CONTENT:                                                                   #
#        defines multiple activation functions for galactic analytical       #
#       potentials.                                                          #
#                                                                            #
# AIM:                                                                       #
#      turn on the potentials gradually so as to avoid artifacts,            #
#   and investigate how different the impacts of them have on the formation  #
#   of spirals.                                                              #
#                                                                            #
# DEPENDENCIES: phycon                                                       #
#                                                                            #
#----------------------------------------------------------------------------#

from src.phycon import *
import math as m

# ~~~~~~~~~~~~~~~~~AF OF HALO~~~~~~~~~~~~~~~~~~~~~
def haloGrow(ti, duration):
    """
    for a halo growth/activation function in reference:
    __________
    calculate the strength of the potential at given time ti with an activation period
    of given duration[Gyr].

    Argument:
    ti: start moment[Gyr] of the potential
    duration: time it takes for the potential to be fully turned on

    Returns:
    potential magnitude at time ti[Gyr]
    """
    ti = ti * gyrTos
    duration = duration * gyrTos
    
    if ti <= duration:
        td = 2.* ti / duration - 1.
        mag = 1. * (3. / 16. * td**5 - 5. / 8. * td**3 + 15. / 16. * td + 1./2.)
    else:
        mag = 1.

    return mag
    
    
# ~~~~~~~~~~~~~~AF OF COX&GOMEZ AMRS~~~~~~~~~~~~~~~
def armEvol(ti, tpeakf = TPEAK / gyrTos, sigma = 1):
    """
    for a arm growth/activation function in reference:
    __________
    calculate the strength of the potential at given time ti with an activation period
    of given duration[Gyr].

    Argument:
    ti: start moment of the potential
    tpeakf: time at which has the peak strength
    sigma: parameter in activiation function suggested by the reference paper

    Returns:
    potential magnitude at time ti[Gyr]
    """
    tpeakf = tpeakf * gyrTos
    sigma = sigma * gyrTos
    ti = ti * gyrTos
    mag = np.exp( - (ti - tpeakf)**2 / (2. * sigma**2) )

    return mag
    
    
# ~~~~~~~~~~~~~~~~AF OF BAR~~~~~~~~~~~~~~~~~~~~~
def barGrowth(ti):
    """
    for a bar growth/activation function in reference:
    __________
    calculate the strength of the potential at given time ti.

    Argument:
    ti: start moment of the potential

    Returns:
    potential magnitude at time ti[Gyr]
    """
    taub = 2 * m.pi / PHIBAR
    t0 = 4 * taub   #time used to reach the half maximum of the total bar strength
    ti = ti * gyrTos
    mag = 1/2 + 1/2 * m.tanh((ti- t0)/taub)

    return mag
    

# ~~~~~~~~~~~~AF OF TIDAL COMPONENT~~~~~~~~~~~~~
def tidalTerm(ti):
    """
    for a tidal component growth/activation function in reference:
    __________
    calculate the strength of the potential at given time ti.

    Argument:
    ti: start moment of the potential
    duration: time it takes for the potential to be fully turned on

    Returns:
    potential magnitude at time ti
    """
    ti = ti * gyrTos
    taut = 1 / PHIBAR_T
    t0 = 2.5 * taut
    mag = m.exp(-(ti - 2 * t0)**2 / (2 * (taut)**2))

    return mag
