# -----------------IMPORT PACAGES--------------------#

from potentials import *
from phycon import *
from grid import *
from spiralFit import *
from coCurve import *
from growth import *

#rotationCurve(addbulge, "black", "Rotation Curve", '+')
os.chdir('/Users/veronicachang/Desktop/pyThesis')
#getDenFileOf("e7bbulge_00000.txt")
#plotPot(cgpot, "potential", "CG Contour Plot", 0.0002)
#---------------get the appendix_pyplot figure in Appendix.-----------
# denAlongR("den/den_bbulge_00050.txt", "#A569BD","default")
# denAlongR("den/den_ebbulge_00000.txt", "#C70039", "s = 3 kpc")
# denAlongR("den/den_e7bbulge_00000.txt", "#7D3C98", "s = 7 kpc")


#---------------get the rotation curve plot.-------------
# rotationCurve(addbulge, "#4A235A", "Total", '-', "mini4rc.txt", 1)
# rotationCurve(halo,"#C70039","Halo",'--')
# rotationCurve(disc, "#3498DB", "Disc", '--')
# rotationCurve(bulge, "#1ABC9C", "Bulge", '--')


#---------------get the contour plot-------------
#comparePot(subPlotStat(haloPot, "potential", 0.), subPlotStat(haloPot, "potential", 0.5),LABEL1="t = 0 Gyr", LABEL2 = "t = 1 Gyr", TITLE="Halo Contour Plot", VMIN = 5*10**10, VMAX = 8.5*10**10)

#comparePot(subPlotStat(cgpot, "potential", 3), subPlotStat(cgpot, "potential", 3.0005),LABEL1="t = 3 Gyr", LABEL2 = "t = 3.0005 Gyr", TITLE="C&G Arms Contour Plot", VMIN = -8 * 10 ** 8, VMAX = 8 * 10 ** 8)
#comparePot(subPlotStat(cgpot, "potential", 2.998), subPlotStat(cgpot, "potential", 2.9995),LABEL1="t = 2.998 Myr", LABEL2 = "t = 2.9985 Myr", TITLE="C&G Arms Contour Plot", VMIN = -8 * 10 ** 8, VMAX = 8 * 10 ** 8)
#comparePot(subPlotStat(testbar, "potential", 0.00000001), subPlotStat(testbar, "potential", 0.08),LABEL1="t = 0 Myr", LABEL2 = "t = 80 Myr", TITLE="Bar Contour Plot", VMIN = -6 * 10 ** 7, VMAX = 6 * 10 ** 7)
