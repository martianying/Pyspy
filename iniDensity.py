import numpy as np
import matplotlib.pyplot as plt
import math
import random 


Rmax = 15

def PDF( r, scale):
    nermerator = np.exp( -r / scale)
    denomenator =  scale - scale * np.exp( - Rmax / scale)#- scale * (np.exp( - Rmax/ scale) - 1)
    return nermerator / denomenator

def CDF( r, scale):
    nermerator = scale - scale * np.exp( - r / scale)
    denomenator = scale - scale * np.exp( - Rmax / scale)
    return nermerator / denomenator

def revCDF( pi, scale):
    part = np.log( (1 - pi) + pi * np.exp(- Rmax / scale))
    return - scale * part
    

#from matplotlib.pyplot import *
ri = np.linspace(0, Rmax, 100)
rii = np.linspace(0, 1, 10000)

rand_basket = np.array([random.random() for i in range(len(rii))])
sorted_random = np.sort(rand_basket)
samples = revCDF( sorted_random, 3)

plt.plot(ri, PDF(ri, 3))
plt.plot(ri, CDF(ri, 3))
plt.hist(samples, bins=Rmax, color='green', density = 1)
#ax = subplot(111)
#ax.set_yscale('log')
#ax.set_xscale('log')
plt.show()