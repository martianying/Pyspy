import numpy as np
import math as m
import matplotlib.pyplot as plt
from matplotlib.pyplot import figure
from numpy import loadtxt
import seaborn as sns
import os
from collections import defaultdict
import shutil

HOME_LOCATION = "/Users/veronicachang/Desktop/pyThesis/datasets"
os.chdir(HOME_LOCATION)
#------To grid the disc
RM = 18                                                  # Radius of the disc 
VEL_UNIT = 10**-5;
RESO_R = 80      
DIS_R = np.linspace(0, RM, RESO_R + 1)                       # Divide R into pieces

RESO_THETA = 60                                          # Divide θ into pieces        
DIS_THETA = np.linspace(0, m.pi/2, RESO_THETA + 1)

PARNUM = 10**6
PARTICLE_MASS = 10**10 / PARNUM 

PAR_LABEL = 0
PAR_V = 1
V_RAD = 2
V_AZI = 3
NUMP = 4

def sectorArea(Rpos):
    width = RM / RESO_R
    length_up = DIS_THETA[1] * DIS_R[Rpos +1] # dθ * R
    length_down = DIS_THETA[1] * DIS_R[Rpos]
    result = 1/ 2 * (length_up + length_down) * width
    return result

def getRloc(xpos, ypos):
    radius = m.sqrt(pow(xpos,2) + pow(ypos, 2))
    Rloc = m.ceil(radius * RESO_R / RM)
    return Rloc

def getVelocity(xVel, yVel):
    Vel = m.sqrt(pow(xVel * VEL_UNIT, 2) + pow(yVel * VEL_UNIT, 2))
    return Vel

def getThetaLoc(xpos, ypos):
    radius = m.sqrt(pow(xpos, 2) + pow(ypos, 2))
    sinValue = abs(ypos / radius)
    cosValue = abs(xpos / radius)

    
    if ((xpos > 0) and (ypos > 0)):
        thetaN = m.ceil(120 * m.asin(sinValue) / m.pi)
    elif ((xpos < 0) and (ypos > 0)) :
        thetaN = m.ceil(120 * m.acos(sinValue) / m.pi) + RESO_THETA
    elif ((xpos < 0) and (ypos < 0)):
        thetaN = m.ceil(120 * m.asin(sinValue) / m.pi) + 2 * RESO_THETA
    else :
        thetaN = m.ceil(120 * m.acos(sinValue) / m.pi) + 3 * RESO_THETA
    
    return thetaN

def othoVelocity( xpos,  ypos,  xVel,  yVel):
    velocityVec = []
    radius = m.sqrt(pow(xpos, 2) + pow(ypos, 2));
    sinValue = ypos / radius;
    cosValue = xpos / radius;
    
    radialV = xVel * VEL_UNIT * cosValue + yVel * VEL_UNIT * sinValue
    velocityVec += [radialV]
    
    azimuthalV =  xVel * VEL_UNIT * sinValue - yVel * VEL_UNIT * cosValue
    velocityVec += [azimuthalV]
    
    return velocityVec


def getPostFile(fileNameIn, fileNameOut):
    fileIn = open(fileNameIn + ".txt", 'r')
    print("---------- Loading ---------")
    fileOut = open(fileNameOut, 'w') 
    while True: 
  
        # Get next line from file 
        line = fileIn.readline()
        if not line: 
            break

        if not line.startswith("#"):
            line = [float(i) for i in line.split()]
            pos1 = str(getRloc(line[0], line[1]))
            pos2 = str(getThetaLoc(line[0], line[1]))
            pos = pos1 + ',' + pos2
    
            vel = getVelocity(line[6], line[7])
            xyvelocity = othoVelocity(line[0], line[1], line[6], line[7])
    
            fileOut.write("%s\t %f\t %f\t %f\t %i\n" % (pos, vel, xyvelocity[0], xyvelocity[1], 1))

    fileIn.close()
    fileOut.close()
    return print("Go check datasets/" + fileNameOut)

def getDenVels(FILE_NAME, paraType):
    if paraType == "DEN":
        CHOICE = NUMP
    if paraType == "VAZI":
        CHOICE = V_AZI
    if paraType == "VRAD":
        CHOICE = V_RAD
    FILE_EXT = paraType
    
    
    g = open(FILE_NAME, 'r')
    lines = g.readlines()
    items= {}

    for line in lines:
        if not line.startswith("#"):
            sline = line.split()
            if sline[PAR_LABEL] in items:
                items[sline[PAR_LABEL]] = float(items.get(sline[PAR_LABEL])) + float(sline[CHOICE])
            else:
                items[sline[PAR_LABEL]] = float(sline[CHOICE])
# Obtained the "items" dictionary
# items: merge values with same key

    newlist = []
    for key,value in items.items():
        keySplited = key.split(",")
        newlist = newlist + [[int(keySplited[0]), int(keySplited[1]), value]]
# Obtaiend the "newlist" list
# newlist: convert the "items" into three dimensional list sat prepare for sorting

    newlist = sorted(newlist, key=lambda x: x[0])
# Obtained the sorted newlist according to its first element

    d = defaultdict(list)
    for k, v1, v2 in newlist:
        d[k].append([v1, v2])
# Subgroupped the newlist, gathered values at same radius
# d is a dictionary

    thetaLocs = RESO_THETA * 4
    allThetaLoc = [i for i in range(thetaLocs)]
    RLocs = RESO_R
    allRloc = [i for i in range(RLocs)]

    dKeys = list(d.keys()) ## avoid overwrite

    for dKey in dKeys:
    #d[dKey] = sorted(d[dKey], key = lambda x: x[0])
        d[dKey] = {i:j for i,j in sorted(d[dKey], key = lambda x: x[0])}
        newSubd = d[dKey]
        absentThetaKeys = set(allThetaLoc).difference(list(newSubd.keys()))
        for absentKey in absentThetaKeys:
            newSubd[absentKey] = 0

    absentRKeys = set(allRloc).difference(dKeys)
    for absentKey in absentRKeys:
        d[absentKey] = {i:0 for i in range(thetaLocs)}

    testFile = open(FILE_NAME[:-4] + FILE_EXT + ".txt", "w") #!!!!!!!!!!
    startKey = 0
    
    while(startKey < RESO_R):
        if CHOICE == NUMP:
            sectorA = sectorArea(startKey) 
        else:
            sectorA = 1
        for subKey in allThetaLoc:
            VALUE = d[startKey][subKey] / sectorA
            testFile.write("%f\t" % (VALUE))
        testFile.write("\n")
        startKey += 1
    testFile.close()
    directory = FILE_NAME.split('/')[0]
    return print(paraType + "file Done. \n Go pyThesis/datasets/" + directory +  "/ check ~" + FILE_EXT + ".txt")
    
def mkdirAndmvFile(FILE_NAME):
    newDirForFile = HOME_LOCATION + "/" + FILE_NAME
    if os.path.exists(newDirForFile):
        message = "directory already exists."
    else: 
        os.mkdir(newDirForFile)
        os.rename(HOME_LOCATION + '/' + FILE_NAME + ".txt", newDirForFile + '/' + FILE_NAME + ".txt")
        message = "ready for further processing."
    return print(message)
    
def desk():
    
    fileNameIn = input("Your input origial Splash file is: ")
    print("~~~~~~~~PROCESSING FOR "+ fileNameIn +".txt STARTS~~~~~~~~")
    mkdirAndmvFile(fileNameIn)
    prefixOfFiles = fileNameIn + "/Post" + fileNameIn
    
    fileNameOut = prefixOfFiles + ".txt"
    fileNameOutDen = prefixOfFiles + "DEN.txt"
    fileNameOutVAZI = prefixOfFiles + "VAZI.txt"
    fileNameOutVRAD = prefixOfFiles + "VRAD.txt"  
    
    
    print("--- All files will be stored in datasets/" + fileNameIn + " ---")
    
    getPostFile(fileNameIn, fileNameOut)

    responce = input("Wanna go ahead? [y / n]")
    print("---------- Loading ---------")
    
    
    if responce == 'y':
        getDenVels(fileNameOut, "DEN")
        print("not yet :P")
        getDenVels(fileNameOut, "VAZI")
        print("not yet :P")
        getDenVels(fileNameOut, "VRAD")
    else:
        print("Post file writing done without furture actions.")
                     
    return print("-------------- FINISHED.---------------")


