#---------------------POTENTIALS-------------------#
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.pyplot import figure
import matplotlib.font_manager
from IPython.core.display import HTML
from matplotlib import rcParams
from matplotlib import ticker
from mpl_toolkits.axes_grid1 import ImageGrid
from mpl_toolkits.axes_grid1 import make_axes_locatable

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
	const = MHALO * Gcode * HALOCI

	xpi = xi * np.cos( OMEH * ti) - yi * np.sin( OMEH * ti)
	ypi = xi * np.sin( OMEH * ti) + yi * np.cos( OMEH * ti)

	rratio = np.sqrt(( ypi**2 + ( xpi**2 + zi**2) * SA**2) / (RHALO**2 * SA**2))
	atanterm = np.arctan( rratio)

	pot = + const * (np.log( rratio) + atanterm / rratio + 0.5 * np.log(( 1. + rratio**2) / ( rratio**2))) / RHALO


	fratio = + const / RHALO * (rratio - np.arctan(rratio)) / ( rratio**3 * RHALO**2)

	fextxpi = - xpi * fratio
	fextypi = - ypi * fratio / SA**2

	fxi = fextxpi * np.cos( -OMEH * ti) - fextypi * np.sin( -OMEH * ti)
	fyi = fextxpi * np.sin( -OMEH * ti) + fextypi * np.cos( -OMEH * ti)
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
	expTerm = np.exp((RARM0 - np.sqrt( xi**2 + yi**2)) / RARMS)

	cirVel = addbulge( ri) * kmTom
	windTerm = cirVel * (ti - TPEAK) / ri

	xpi = xi * np.cos( windTerm) - yi * np.sin( windTerm)
	ypi = yi * np.cos( windTerm) + xi * np.sin( windTerm)

	n_list = [1, 2, 3]

	sumPart = 0

	for i in range( len( n_list) ):
		Kn = n_list[i] * NARM / (np.sqrt( xpi**2 + ypi**2) * np.sin( PITCHARM))
		Bn = Kn * HARM * (1 + 0.4 * Kn * HARM)
		# no sechTerm here sat simplify the problem
		Dn = (1 + Kn * HARM + 0.3 * (Kn * HARM)**2) / (1 + 0.3 * Kn * HARM)
		gamma = NARM * (np.arctan2( ypi, xpi) - np.log(np.sqrt( xpi**2 + ypi**2) / RARM0) / np.tan(PITCHARM))

		sumPart = sumPart + CnARM[i] / ( Kn * Dn) * np.cos(( i + 1) * gamma)

	pot = As * CGARM * expTerm * sumPart

	return pot



#----bar
def testbar(xi, yi, zi, ti):
	xi = xi * kpcTom
	yi = yi * kpcTom
	zi = zi * kpcTom
	ti = ti * gyrTos


	xpi = xi * np.cos( ti * PHIBAR) - yi * np.sin( ti * PHIBAR)
	ypi = xi * np.sin( ti * PHIBAR) + yi * np.cos( ti * PHIBAR)

	arcterm = 2 * np.arctan2( ypi, xpi )
	five = ( LMABAR * BETA2)**5
	r = np.sqrt( xi**2 + yi**2 + zi**2)
	denofive = (1. + r**5 / five)
   
	const = - barGrowth(ti) * Gcode * MBAR * ALPHA2 / LMABAR**3
   
	d2 = xi**2 + yi**2
	pot = const * d2 * np.cos(arcterm) / denofive
   
	fextxpi = const * (5. * xpi * d2 * r**3 * np.cos(arcterm) / (five * denofive **2) - (2. * xpi * np.cos(arcterm) + 2. * ypi * np.sin(arcterm))/denofive)

	fextypi = const * (5. * ypi * d2 * r**3 * np.cos(arcterm) / (five * denofive**2) - (2. * ypi * np.cos(arcterm) - 2. * xpi * np.sin(arcterm))/denofive)

	fxi = + (fextxpi * np.cos( -PHIBAR * ti) - fextypi * np.sin( -PHIBAR * ti))
	fyi = + (fextxpi * np.sin( -PHIBAR * ti) + fextypi * np.cos( -PHIBAR * ti))
	fzi = const * 5. * zi * d2 * r**3 * np.cos(arcterm) / (five * denofive**2)

	return pot, fxi, fyi


#-----------------PLOT FUNCTION-----------------#


def subPlotStat(func, porf, time):
	x, y = np.linspace(-15, 15, 600), np.linspace(-15, 15, 600)
	X, Y = np.meshgrid(x, y)
	Z = np.zeros((600, 600))

	if porf == "potential":
		if func == haloPot:
			T = func(X, Y, Z, time, 0.4, 1)[0]
		elif func == cgpot:
			T = func(X, Y, Z, time)
		else:
			T = func(X, Y, Z, time)[0]
	elif port == "force":
		if func == haloPot:
			T = func(X, Y, Z, time, 0.5, 3)[1]
		else:
			T = func(X, Y, Z, time)[1]

	return X, Y, T




def comparePot(sub_plot_stat_1, sub_plot_stat_2, TITLE = "title", LABEL1="time1", LABEL2="time2", VMIN=-8 * 10 ** 8,
			   VMAX=8 * 10 ** 8):
	fig = plt.figure(figsize=(8, 4.5))

	grid = ImageGrid(fig, 111,  # as in plt.subplot(111)
					 nrows_ncols=(1, 2),
					 axes_pad=0.052,  # gap between subplots
					 share_all=True)

	fig.suptitle(TITLE, fontsize=12.5, x=0.5, y=0.92)

	ax1, ax2 = grid[0], grid[1]

	im1 = ax1.imshow(sub_plot_stat_1[2], extent=[-15, 15, -15, 15], cmap='RdGy', vmin=VMIN, vmax=VMAX)
	im2 = ax2.imshow(sub_plot_stat_2[2], extent=[-15, 15, -15, 15], cmap='RdGy', vmin=VMIN, vmax=VMAX)

	for ax in grid:
		ax.set_xticks([-10, -5, 0, 5, 10])
		ax.set_yticks([-10, -5, 0, 5, 10])
		ax.set_xlabel('X [kpc]')
		ax.set_ylabel('Y [kpc]')
		ax.tick_params(direction="in")

	ax1.text(-14.5, 13.64, LABEL1, bbox={'facecolor': 'white', 'pad': 2.7})
	ax2.text(-14.5, 13.64, LABEL2, bbox={'facecolor': 'white', 'pad': 2.7})

	cax = plt.axes([0.908, 0.154, 0.03, 0.682])  # adjust the axes positon (left, bottom, width, height)
	cbar = fig.colorbar(im2, cax=cax)  #
	cbar.ax.tick_params(labelsize=8, grid_alpha=0.5, direction='in')
	cbar.set_label(label=r'$cm^{2}/s^{2}$', size='medium')

	return plt.savefig('images/'+ TITLE + '.png', dpi=300)


