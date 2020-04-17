#---------------------POTENTIALS-------------------#
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.pyplot import figure
import matplotlib.font_manager
from IPython.core.display import HTML
from matplotlib import rcParams
from matplotlib import ticker

from phycon import *
from roCurve import *
from growth import *


#----halo
def haloPot(xi, yi, zi, ti, q, span):
	xi = xi * kpcTom
	yi = yi * kpcTom
	zi = zi * kpcTom
	ti = ti * gyrTos

	SA = 1. - haloGrow(ti / gyrTos, span) * q
	const = Mhalo * Gcode * HaloCi

	xpi = xi * np.cos( omeH * ti) - yi * np.sin( omeH * ti)
	ypi = xi * np.sin( omeH * ti) + yi * np.cos( omeH * ti)

	rratio = np.sqrt(( ypi**2 + ( xpi**2 + zi**2) * SA**2) / (rhalo**2 * SA**2))
	atanterm = np.arctan( rratio)

	pot = + const * (np.log( rratio) + atanterm / rratio + 0.5 * np.log(( 1. + rratio**2) / ( rratio**2))) / rhalo


	fratio = + const / rhalo * (rratio - np.arctan(rratio)) / ( rratio**3 * rhalo**2)

	fextxpi = - xpi * fratio
	fextypi = - ypi * fratio / SA**2

	fxi = fextxpi * np.cos( -omeH * ti) - fextypi * np.sin( -omeH * ti)
	fyi = fextxpi * np.sin( -omeH * ti) + fextypi * np.cos( -omeH * ti)
	fzi = - zi * fratio

	return pot, fxi, fyi



#----cgArm
def cgpot(xi, yi, zi, ti):
	xi = xi * kpcTom
	yi = yi * kpcTom
	zi = zi * kpcTom
	ti = ti * gyrTos
	ri = np.sqrt(xi**2 + yi**2 + zi**2)

	As = armEvol(ti / gyrTos)
	expTerm = np.exp((r0 - np.sqrt( xi**2 + yi**2)) / Rs)

	cirVel = addbulge( ri) * kmTom
	windTerm = cirVel * (ti - tpeak) / ri

	xpi = xi * np.cos( windTerm) - yi * np.sin( windTerm)
	ypi = yi * np.cos( windTerm) + xi * np.sin( windTerm)

	n_list = [1, 2, 3]

	sumPart = 0

	for i in range( len( n_list) ):
		Kn = n_list[i] * Nm / (np.sqrt( xpi**2 + ypi**2) * np.sin( pA))
		Bn = Kn * H * (1 + 0.4 * Kn * H)
		# no sechTerm here sat simplify the problem
		Dn = (1 + Kn * H + 0.3 * (Kn * H)**2) / (1 + 0.3 * Kn * H)
		gamma = Nm * (np.arctan2( ypi, xpi)
						 - np.log(np.sqrt( xpi**2 + ypi**2) / r0) / np.tan(pA))

		sumPart = sumPart + Cn[i] / ( Kn * Dn) * np.cos(( i + 1) * gamma)

	pot = As * constCG * expTerm * sumPart

	return pot

#----bar
def testbar(xi, yi, zi, ti):
	xi = xi * kpcTom
	yi = yi * kpcTom
	zi = zi * kpcTom


	xpi = xi * np.cos( ti * phibar) - yi * np.sin( ti * phibar)
	ypi = xi * np.sin( ti * phibar) + yi * np.cos( ti * phibar)

	arcterm = 2 * np.arctan2( ypi, xpi )
	five = ( LMabar * beta2)**5
	r = np.sqrt( xi**2 + yi**2 + zi**2)
	denofive = (1. + r**5 / five)
   
	const = - barGrowth(ti) * Gcode * Mbar * alpha2 / LMabar**3
   
	d2 = xi**2 + yi**2
	pot = const * d2 * np.cos(arcterm) / denofive
   
	fextxpi = const * (5. * xpi * d2 * r**3 * np.cos(arcterm) / (five * denofive **2) - (2. * xpi * np.cos(arcterm) + 2. * ypi * np.sin(arcterm))/denofive)

	fextypi = const * (5. * ypi * d2 * r**3 * np.cos(arcterm) / (five * denofive**2) - (2. * ypi * np.cos(arcterm) - 2. * xpi * np.sin(arcterm))/denofive)

	fxi = + (fextxpi * np.cos( -phibar * ti) - fextypi * np.sin( -phibar * ti))
	fyi = + (fextxpi * np.sin( -phibar * ti) + fextypi * np.cos( -phibar * ti))
	fzi = const * 5. * zi * d2 * r**3 * np.cos(arcterm) / (five * denofive**2)

	return pot, fxi, fyi


#-----------------PLOT FUNCTION-----------------#
#x, y = np.linspace(-10, 10, 100), np.linspace(-10, 10, 100)
#X, Y = np.meshgrid(x, y)
#Z = np.zeros((100,100))


def plotPot(func, porf, title):

	x, y = np.linspace(-10, 10, 100), np.linspace(-10, 10, 100)
	X, Y = np.meshgrid(x, y)
	Z = np.zeros((100,100))
	
	if porf == "contour":
		T = func(X, Y, Z, 2.35, 0.5, 3)[0]
	else:
		T = func(X, Y, Z, 2.35, 0.5, 3)[1]

	figure(num=None, figsize=(2, 1.8), dpi=300, facecolor='w', edgecolor='k')
	font = {'family' : 'DejaVu Sans', 'size': 8}
	mpl.rc('font', **font)              # pass in the font dict as kwargs
	rcParams['axes.titlepad'] = 12

	ax = plt.axes()
	ax.tick_params( direction="in")

	plt.contourf(X, Y, T, 20, cmap='RdGy')

	#plt.clim(-2*10**8,2*10**8)
	plt.title(title + " Contour Plot")

	plt.xlabel("x / kpc")
	plt.xticks(np.arange(-10, 15, step = 5),fontsize = 5,fontweight = 'normal')

	plt.ylabel("y / kpc")
	plt.yticks(np.arange(-10, 15, step = 5),fontsize = 5, fontweight='normal')


	tick_locator = ticker.MaxNLocator(nbins = 5)
	cbar = plt.colorbar( extend = 'neither', spacing = 'proportional',
				orientation = 'vertical', shrink = 1)
	cbar.locator = tick_locator
	cbar.update_ticks()
	cbar.ax.tick_params(labelsize = 5, grid_alpha = 0.5, direction = 'in')
	cbar.set_label(label = r'$cm^{2}/s^{2}$', size = 'small')

	return cbar

# #plotPot(haloPot, "force","halo")
# #plt.show()  

