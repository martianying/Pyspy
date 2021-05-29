##################################   PLAYGOURND   ######################################
# you could play with the code below to generate images I have put into my thesis
# Unfortunately, I could not provide any data from the research at this moment since it might
# will be published as a paper in the future.

import os
from src.roCurve import *
from src.potentials import *

# set working directory
os.chdir('/Users/.../Desktop/...')

if __name__ == '__main__':
    # Plot the rotation curve
    rotationCurve(addbulge, "black", "Rotation Curve", '+')
