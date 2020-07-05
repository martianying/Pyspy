import numpy as np
import math as m
import matplotlib.pyplot as plt
from matplotlib.pyplot import figure
from numpy import loadtxt
import seaborn as sns
import os
from collections import defaultdict
import shutil
import matplotlib.cm as cm
from scipy.optimize import curve_fit

HOME_LOCATION = "/Users/veronicachang/Desktop/pyThesis/datasets"
os.chdir(HOME_LOCATION)
# ------To grid the disc
RM = 15  # Radius of the disc
VEL_UNIT = 10 ** -5;
RESO_R = 80
DIS_R = np.linspace(0, RM, RESO_R + 1)  # Divide R into pieces

RESO_THETA = 100  # Divide θ into pieces
DIS_THETA = np.linspace(0, m.pi / 2, RESO_THETA + 1)

PARNUM = 10 ** 6
PARTICLE_MASS = 10 ** 10 / PARNUM

PAR_LABEL = 0
PAR_V = 1
V_RAD = 2
V_AZI = 3
NUMP = 4

def getNumOfLineR(rKpc):
    return m.ceil(rKpc * RESO_R / RM)

low = getNumOfLineR(2)
high = getNumOfLineR(7)

def sectorArea(Rpos):
    width = RM / RESO_R
    length_up = DIS_THETA[1] * DIS_R[Rpos + 1]  # dθ * R
    length_down = DIS_THETA[1] * DIS_R[Rpos]
    result = 1 / 2 * (length_up + length_down) * width
    return result * 10**6


def getRloc(xpos, ypos):
    radius = m.sqrt(pow(xpos, 2) + pow(ypos, 2))
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
        thetaN = m.ceil(2 * RESO_THETA * m.asin(sinValue) / m.pi)
    elif ((xpos < 0) and (ypos > 0)):
        thetaN = m.ceil(2 * RESO_THETA * m.acos(sinValue) / m.pi) + RESO_THETA
    elif ((xpos < 0) and (ypos < 0)):
        thetaN = m.ceil(2 * RESO_THETA * m.asin(sinValue) / m.pi) + 2 * RESO_THETA
    else:
        thetaN = m.ceil(2 * RESO_THETA * m.acos(sinValue) / m.pi) + 3 * RESO_THETA

    return thetaN


def othoVelocity(xpos, ypos, xVel, yVel):
    velocityVec = []
    radius = m.sqrt(pow(xpos, 2) + pow(ypos, 2));
    sinValue = ypos / radius;
    cosValue = xpos / radius;

    radialV = xVel * VEL_UNIT * cosValue + yVel * VEL_UNIT * sinValue
    velocityVec += [radialV]

    azimuthalV = xVel * VEL_UNIT * sinValue - yVel * VEL_UNIT * cosValue
    velocityVec += [azimuthalV]

    return velocityVec


def getPostFile(fileNameIn, fileNameOut):
    fileIn = open(fileNameIn + ".txt", 'r')
    i = 0
    while (i < 14):
        fileIn.readline()
        i = i + 1

    print("---------- Loading ---------")
    fileOut = open(fileNameOut, 'w')
    fileOut.write("#This is the processed ascii file for " + fileNameIn + ".txt \n")
    fileOut.write("#{Rloc,thetaLoc} \t (velocity km/s) \t (radial velocity) \t (azimuthal velocity) \t (#particle) \n")

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
    if paraType == "NUM":
        CHOICE = NUMP
    if paraType == "VAZI":
        CHOICE = V_AZI
    if paraType == "VRAD":
        CHOICE = V_RAD
    FILE_EXT = paraType

    g = open(FILE_NAME, 'r')
    lines = g.readlines()
    items = {}

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
    for key, value in items.items():
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

    dKeys = list(d.keys())  ## avoid overwrite

    for dKey in dKeys:
        # d[dKey] = sorted(d[dKey], key = lambda x: x[0])
        d[dKey] = {i: j for i, j in sorted(d[dKey], key=lambda x: x[0])}
        newSubd = d[dKey]
        absentThetaKeys = set(allThetaLoc).difference(list(newSubd.keys()))
        for absentKey in absentThetaKeys:
            newSubd[absentKey] = 0

    absentRKeys = set(allRloc).difference(dKeys)
    for absentKey in absentRKeys:
        d[absentKey] = {i: 0 for i in range(thetaLocs)}

    testFile = open(FILE_NAME[:-4] + FILE_EXT + ".txt", "w")  # !!!!!!!!!!
    startKey = 0

    while (startKey < RESO_R):
        for subKey in allThetaLoc:
            VALUE = d[startKey][subKey]
            testFile.write("%f\t" % (VALUE))
        testFile.write("\n")
        startKey += 1
    testFile.close()
    directory = FILE_NAME.split('/')[0]
    return print(paraType + ".txt done.")


def mkdirAndmvFile(FILE_NAME):
    oldDirForFile = HOME_LOCATION + '/'
    oldFileName = HOME_LOCATION + '/' + FILE_NAME + ".txt"

    newDirForFile = HOME_LOCATION + "/" + FILE_NAME
    newFileName = HOME_LOCATION + "/" + FILE_NAME + "/" + FILE_NAME + ".txt"

    if os.path.exists(newDirForFile):
        message = "directory already exists."
    else:
        os.mkdir(newDirForFile)
    if os.path.isfile(oldFileName):
        shutil.move(oldFileName, newFileName)
    message = "ready for further processing."
    return print(message)


def rmmvPostFiles(DIRE):
    toRemoveVAZI = HOME_LOCATION + "/" + DIRE + "/Post" + DIRE + "VAZI.txt"
    toRemoveVRAD = HOME_LOCATION + "/" + DIRE + "/Post" + DIRE + "VRAD.txt"
    toRenameVAZIaved = HOME_LOCATION + "/" + DIRE + "/Post" + DIRE + "VAZIaved.txt"
    toRemoveVRADaved = HOME_LOCATION + "/" + DIRE + "/Post" + DIRE + "VRADaved.txt"

    os.remove(toRemoveVAZI)
    os.remove(toRemoveVRAD)

    shutil.move(toRenameVAZIaved, toRenameVAZIaved[:-8] + ".txt")
    shutil.move(toRemoveVRADaved, toRemoveVRADaved[:-8] + ".txt")

    return print("finish clean-up")


def averageVFile(PostVFileName, PostNumFileName):
    numfile = open(PostNumFileName, 'r')
    vfile = open(PostVFileName, 'r')
    aVfile = open(PostVFileName[:-4] + "aved.txt", 'w')
    i = 0
    while (i < RESO_R):
        lineNum = np.array(numfile.readline().split()).astype(float) + 0.0001
        lineV = np.array(vfile.readline().split()).astype(float)
        averagedV = lineV / lineNum
        if PostVFileName[-8:-4] == "VAZI":
            constrast = np.average(averagedV)
            putin = averagedV - constrast
        else:
            putin = averagedV
        for velInGrid in putin:
            aVfile.write("%f\t" % (velInGrid))
        aVfile.write("\n")
        i += 1

    numfile.close()
    vfile.close()
    aVfile.close()
    return print("--> Obtained aved file")


def getDENfileFromNUM(PostNumFileName):
    numfile = open(PostNumFileName, 'r')
    denfile = open(PostNumFileName[:-7] + "DEN.txt", 'w')
    i = 0
    while (i < RESO_R):
        lineNum = np.array(numfile.readline().split()).astype(float)
        sectorA = sectorArea(i)
        for numP in lineNum:
            gridDen = numP * PARTICLE_MASS/ sectorA
            denfile.write("%f\t" % (gridDen))
        denfile.write("\n")
        i += 1

    numfile.close()
    denfile.close()
    return print("--> Obtained DEN file")


def desk():
    fileNameIn = input("Your input origial Splash file is: ")
    print("~~~~~~~~PROCESSING FOR " + fileNameIn + ".txt STARTS~~~~~~~~")
    mkdirAndmvFile(fileNameIn)
    prefixOfOutFiles = fileNameIn + "/Post" + fileNameIn

    fileNameOut = prefixOfOutFiles + ".txt"
    fileNameOutNUM = prefixOfOutFiles + "NUM.txt"
    fileNameOutVAZI = prefixOfOutFiles + "VAZI.txt"
    fileNameOutVRAD = prefixOfOutFiles + "VRAD.txt"

    print("--- All files will be stored in datasets/" + fileNameIn + "/ ---")

    getPostFile(fileNameIn + '/' + fileNameIn, fileNameOut)

    responce = input("Wanna go ahead? [y / n]")
    print("---------- Loading ---------")

    if responce == 'y':
        getDenVels(fileNameOut, "NUM")
        print("not yet :P")
        getDenVels(fileNameOut, "VAZI")
        averageVFile(fileNameOutVAZI, fileNameOutNUM)
        print("not yet :P")
        getDenVels(fileNameOut, "VRAD")
        averageVFile(fileNameOutVRAD, fileNameOutNUM)
        print("not yet :P")
        getDENfileFromNUM(fileNameOutNUM)

        rmmvPostFiles(fileNameIn)

    else:
        print("Post file writing done without furture actions.")

    return print("-------------- FINISHED.---------------")


import matplotlib.cm as cm


def plotPostFileOf(FILE_NAME, fileType, TITLE="", xAxis=1, pos=''):
    Plot_FILE_NAME = FILE_NAME + "/Post" + FILE_NAME + fileType + ".txt"
    heatMapData = open(Plot_FILE_NAME, 'r')
    lines = heatMapData.readlines()
    loaded_data = []
    for line in lines:
        sline = line.split()
        loaded_data += [sline]
    dendata = np.array(loaded_data).astype(float)

    if fileType == "VAZI":
        Vmax = 12
        Vmin = - Vmax
        colorBarComment = r'$ \langle v_{\phi} - \overline{v}_{\phi} \rangle \, [km/s]$'
    elif fileType == "DEN":
        Vmax = 1300
        Vmin = 0
        colorBarComment = r'$Σ_{g} \, [\log{(M_{\odot}pc^{-2})}$]'
    else:
        Vmax = 20
        Vmin = -Vmax
        colorBarComment = r'$v_{R} \, [km/s]$'

    heat_map = sns.heatmap(dendata, cmap=cm.coolwarm,
                           cbar_kws={'label': colorBarComment, "shrink": .9}, vmin=Vmin, vmax=Vmax)
    # ax.txt(fontsize=12)
    cbar = heat_map.collections[0].colorbar
    cbar.ax.tick_params(labelsize=5, direction="in")

    if fileType == "DEN":
        cbar.set_ticks([10, 100, 1000])
        # cbar.set_ticks([0, 10, 100, 1000])
        cbar.set_ticklabels(['1','2', '3'])

    if fileType == "VRAD":
        cbar.set_ticks([Vmin, Vmin * 1 / 2, 0, Vmax * 1 / 2, Vmax])

    if pos == "middle":
        cbar.set_ticks([Vmin, Vmin * 1 / 2, 0, Vmax * 1 / 2, Vmax])
        plt.ylabel("R [kpc]")

    plt.xticks([0, RESO_THETA, 2 * RESO_THETA, 3 * RESO_THETA, 4 * RESO_THETA],
               ['0', '$\pi$ / 2', '$\pi$', '$3\pi / 2$', "2$\pi$"], rotation=0, fontsize=7)

    plt.yticks([0, RESO_R * 1 / 3, RESO_R * 2 / 3],
               ['0', '5', '10'], rotation=0, fontsize=7)

    plt.xlabel(r"$\theta$ [rad]")

    if xAxis == 0:
        plt.tick_params(
            axis='x',  # changes apply to the x-axis
            which='both',  # both major and minor ticks are affected
            bottom=True,  # ticks along the bottom edge are off
            top=False,  # ticks along the top edge are off
            labelbottom=False)  # labels along the bottom edge are off
        plt.xlabel("")

    plt.tick_params(direction="in")
    plt.title(TITLE, x=0.5, y=1.05)
    # plt.xlabel(r"$\theta$ [rad]")

    return plt.savefig('../images/' + "TEST" + '.png', dpi=300)


def plotPostFilesPanel():
    dirName = input("Enter oritginial Splash file name: ")
    heatMapName = input("Name the Heatmap: ")
    plt.subplots_adjust(wspace=0, hspace=0.08)

    plt.subplot(311)
    plotPostFileOf(dirName, "DEN", heatMapName, xAxis=0)

    plt.subplot(312)
    plotPostFileOf(dirName, "VAZI", xAxis=0, pos="middle")

    plt.subplot(313)
    plotPostFileOf(dirName, "VRAD", xAxis=1)

    plt.savefig(HOME_LOCATION + '/' + dirName + '/' + dirName + "heatMap.png", dpi=300)

    return print("images grab go datasets/" + dirName + "/ check heatMap.png")


def logPostDenHeatmap(dirName):
    title = input("title of the plot: ")
    import matplotlib.cm as cm

    # cb = plt.colorbar(ticks=[1,5,10,20,50], format=formatter)
    # ax = fig.add_subplot(111)
    plt.tick_params(direction='in', which='both')

    heatMapData = open(dirName + '/Post' + dirName + "DEN.txt", 'r')
    lines = heatMapData.readlines()

    y, x = np.linspace(0, RM, RESO_R), np.linspace(0, 2 * m.pi, 4 * RESO_THETA)
    zValues = []
    for line in lines:
        row = line.split()
        zValues += [row]
    np.array(zValues).astype(float)

    Vmax = 1300.
    plt.contourf(x, y, zValues, 20, cmap=cm.coolwarm, vmin=0., vmax=Vmax)
    plt.yscale('log')
    cb = plt.colorbar()

    cb.ax.tick_params(axis='y', direction='in')
    #cb.set_label(r'$Σ_{g} \, [\log{(M_{\odot}pc^{-2})}$]')
    cb.set_label(r'$Σ_{g} \, [M_{\odot}pc^{-2}$]')
    cb.set_ticks([10, 100, 1000])
    cb.set_ticklabels(['10', '100', '1000'])

    plt.title(title)

    plt.xticks([0, m.pi / 2, m.pi, m.pi * 3 / 2, 2 * m.pi],
               ['0', '$\pi$ / 2', '$\pi$', '$3\pi / 2$', "2$\pi$"], rotation=0)

    plt.ylim(0.9, RM)
    plt.ylabel("log(R) [kpc]")
    plt.xlabel(r"$\theta$ [rad]")
    plt.gca().invert_yaxis()
    return plt.savefig(HOME_LOCATION + '/' + dirName + '/' + dirName + "DENlogged.png", dpi=300)


def getRlineIn(fileName, radiusR):
    g = open(fileName, 'r')
    i = 0
    getLineNum = m.ceil(radiusR / RM * RESO_R)
    while(i < getLineNum):
        rline = g.readline()
        i += 1
    g.close()
    rline = np.array(rline.split()).astype(float)
    return rline


def trigoPlotOf(dirName, R=3):
    xaxis = np.linspace(0, 2 * m.pi, 4 * RESO_THETA)[1:]

    prefix = dirName + '/Post' + dirName
    data1 = getRlineIn(prefix + "DEN.txt", R)[1:]  # [1:] not [:] because of m.ceil
    data2 = getRlineIn(prefix + "VAZI.txt", R)[1:]
    data3 = getRlineIn(prefix + "VRAD.txt", R)[1:]

    fig = plt.figure()
    ax = fig.add_subplot(111)

    plt.yscale('log')
    ax.tick_params(axis='both', which='both', direction='in')
    denPlot = ax.plot(xaxis, data1, '-r', label=r'$Σ_{g}$')

    ax2 = ax.twinx()
    ax2.tick_params(axis='y', direction='in')
    vradPlot = ax2.plot(xaxis, data3, '-', label=r'$v_{R}$')
    vaziPlot = ax2.plot(xaxis, data2, '-', label=r'$\langle v_{\theta} - \overline{v}_{\theta} \rangle$')

    # added these three lines
    lns = denPlot + vradPlot + vaziPlot
    labs = [l.get_label() for l in lns]
    ax.legend(lns, labs, loc=1)

    # ax.grid()
    plt.xticks([0, m.pi / 2, m.pi, 3 / 2 * m.pi, 2 * m.pi], ['0', '$\pi$ / 2', '$\pi$', '$3\pi / 2$', "2$\pi$"],
               rotation=0, fontsize=7)
    ax.set_xlabel(r"$\theta$ [rad]")
    ax.set_ylabel(r'$Σ_{g} \, [(M_{\odot}pc^{-2})$]')
    ax2.set_ylabel("V [km/s]")

    ax.set_ylim(100, 2000)
    ax2.set_ylim(-50, 50)

    return plt.savefig(HOME_LOCATION + '/' + dirName + '/' + dirName + "trigo" + str(R) + "kpc.png", dpi=300)


from scipy.optimize import curve_fit
xdata = np.linspace(0, 2 * m.pi, 4 * RESO_THETA)[1:]


def FourierFunc(x, a0, a1, a2, a3, a4, a5, a6, w1, w2, w3, w4, w5, w6):
    '''
    Fourier serise model for fitting
    '''
    part1 = a0 + a1 * np.cos(x - w1) + a2 * np.cos(2 * (x - w2)) + a3 * np.cos(3 * (x - w3))
    part2 = a4 * np.cos(4 * (x - w4)) + a5 * np.cos(5 * (x - w5)) + a6 * np.cos(6 * (x - w6))
    result = part1 + part2
    return result


def plotFourierFit(dirName, rKpc):
    fileName = dirName + '/Post' + dirName + 'DEN.txt'
    ydata = getRlineIn(fileName, rKpc)[1:]

    popt, pcov = curve_fit(FourierFunc, xdata, ydata)

    paraValues = popt
    paraKeys = ['a0','a1','a2','a3','a4','a5','a6','w1','w2','w3','w4','w5','w6']
    result = dict(zip(paraKeys, paraValues))
    result = paraValues
    plt.plot(xdata, ydata, label='data', linewidth = 4, alpha = 0.5)
    plt.plot(xdata, FourierFunc(xdata, *popt),  label='fit', linewidth = 3, alpha = 50, color = "black")
    return plt.show()


def FourierFunc(x, a0, a2, a4, a6, w2, w4, w6):
    '''
    Fourier serise model for fitting
    '''
    part1 = a0 + a2 * np.cos(2 * (x - w2))
    part2 = a4 * np.cos(4 * (x - w4)) + a6 * np.cos(6 * (x - w6))
    result = part1 + part2
    return result


def getAllFitPhaseof(dirName, fileType, iPhase):
    fileName = dirName + "/Post" + dirName + fileType + ".txt"
    xdata = np.linspace(0, 2 * m.pi, 4 * RESO_THETA)[1:]
    result = []
    i = 0
    while (i < RESO_R):
        radius = (i + 1) * RM / RESO_R
        ydata = getRlineIn(fileName, radius)[1:]

        popt, pcov = curve_fit(FourierFunc, xdata, ydata)

        paraValues = popt
        paraKeys = ['a0', 'a2', 'a4', 'a6', 'w2', 'w4', 'w6']
        paraDict = dict(zip(paraKeys, paraValues))
        result += [paraDict[iPhase]]
        paraDict.clear()
        i += 1

    return np.array(result[low:high])


def adjustiPhase(arriPhase):
    offset = (np.abs(np.min(arriPhase)) // np.pi) * np.pi + np.pi
    adjusted = arriPhase + offset
    yis = adjusted
    i = 0
    while i < len(yis) - 1:
        if yis[i] < yis[i + 1] and (yis[i] + m.pi) > yis[i + 1]:
            i += 1
        else:
            if yis[i] > yis[i + 1]:
                yis[i + 1] += m.pi
            else:
                yis[i + 1] -= m.pi

    result = yis - yis[0]
    return result


def getNewPhases(data):
    peorids = np.array(data) // np.pi
    newPhases = np.array(data) - peorids * np.pi
    return newPhases


def comparePhasesPot():
    dirName = input("directory you want to work on: ")
    print("choose a phase to plot:\n a)w2 b) w4 c) w6 \n")
    choice = input("your choice is: ")
    plotName = input("name the plot: ")

    if choice == 'a':
        iPhase = 'w2'
    if choice == 'b':
        iPhase = 'w4'
    if choice == 'c':
        iPhase = 'w6'

    denDataset = getAllFitPhaseof(dirName, 'DEN', iPhase)
    vAZIDataset = getAllFitPhaseof(dirName, 'VAZI', iPhase)
    vRADDataset = getAllFitPhaseof(dirName, 'VRAD', iPhase)

    rdata = np.linspace(0, RM, RESO_R)[low:high]

    plt.tick_params(axis='both', direction='in')

    plt.plot(rdata, adjustiPhase(denDataset), label=r'$Σ_{g}$', color='r')
    plt.plot(rdata, adjustiPhase(vAZIDataset), label=r'$ \langle v_{\theta} - \overline{v}_{\theta} \rangle $')
    plt.plot(rdata, adjustiPhase(vRADDataset), label=r'$v_{R}$')

    plt.title(plotName)
    plt.ylabel(r'$\varphi \, [rad]$')
    plt.xlabel(r'$R \, [kpc]$')
    plt.legend(loc="upper left")

    return plt.savefig(HOME_LOCATION + '/' + dirName + '/' + dirName + iPhase + ".png", dpi=300)

