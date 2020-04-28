# -----------------IMPORT PACAGES--------------------#
import numpy as np
import math as m
import matplotlib as mpl
import matplotlib.pyplot as plt

from matplotlib.pyplot import figure
import matplotlib.font_manager
from IPython.core.display import HTML

from matplotlib import rcParams
from matplotlib import ticker

from potentials import *
from phycon import *
from grid import *

#rotationCurve(addbulge, "black", "Rotation Curve", '+')
os.chdir('/Users/veronicaplanck/Desktop/appendixA')
getDenFileOf("blank80.txt")
