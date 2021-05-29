#----------------------------------------------------------------------------#
# MODULE: spiralFit                                                          #
#                                                                            #
# CONTENT:                                                                   #
#      Perform the Fourier fitting.                                          #
#                                                                            #
# AIM:                                                                       #
#      with this module, we could plot figure 3.14 in my thesis, which       #
#   address the phase change and relation between density, azimuthal and     #
#   radial velocity of gas particles.                                        #
#                                                                            #
#                                                                            #
# REFERENCE:                                                                 #
# DEPENDENCIES: grid                                                         #
#                                                                            #
#----------------------------------------------------------------------------#
from src.grid import *
from scipy.optimize import curve_fit

#~~~~~~~~~~~~~~~~~~~~~~~~~~~ FOURIER FITTING~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# model 
def FourierFunc(x, a0, a1, a2, a3, a4, a5, a6, w1, w2, w3, w4, w5, w6):
    '''
    HELPER function
    Fourier serise model for fitting

    Arguments:
    x:       data to be fitted
    others:  parameters of the model

    Returns:
    parameters
    '''
    part1 = a0 + a1 * np.cos(x - w1) + a2 * np.cos(2 * (x - w2)) + a3 * np.cos(3 * (x - w3))
    part2 = a4 * np.cos(4 * (x - w4)) + a5 * np.cos(5 * (x - w5)) + a6 * np.cos(6 * (x - w6))
    result = part1 + part2

    return result


def fourierFit(denFileName, Rloc = 40):
    '''
    Fit density series at each radius with:
    'y' for printing parameters dictionary
    'n' for plotting fitting curve

    Arguments:
    denFileName: file stores number density to be Fourier fitted
    Rloc: index of the radius
    (NOTE: adjusted_denfile is one line in that file instead of the whole)

    Returns:
    plot and store the scattered number density and its fourier fit.
    '''
    CHOICE = input("Print Parameters (Y/y) or Plot (N/n): ") or "Y"
    xdata = np.linspace(0, 2 * m.pi, RESO_THETA * 4) 
    ydata = adjustFile(denFileName, 'y')[Rloc]

    # perform the fitting
    popt, pcov = curve_fit(FourierFunc, xdata, ydata)

    if CHOICE == 'y' or CHOICE == 'Y':
        paraValues = np.abs(popt)
        paraKeys = ['a0','a1','a2','a3','a4','a5','a6','w1','w2','w3','w4','w5','w6']
        result = dict(zip(paraKeys, paraValues))
        print("----------------------------")
        print("Fourier fit at R = " + str(round(Rloc * DIS_R[1], 2)) + " kpc is :")
        print(result)
        print("----------------------------")
    
    if CHOICE == 'n' or CHOICE == 'N':
        (figure, axes) = plt.subplots()
        axes.yaxis.set_tick_params(direction='in', which='both')
        axes.xaxis.set_tick_params(direction='in', which='both')

        plt.plot(xdata, ydata, label='data', color = COLOR1, linewidth = 4, alpha = 0.5)
        plt.plot(xdata, FourierFunc(xdata, *popt),  label='fit', linewidth = 3, alpha = 50, color = "black")
        plt.xlabel(r'$\theta$ [rad]')
        plt.ylabel(r'Contrast Number Density')
        plt.legend(loc = "lower right")
        plt.title("Fourier Fit R = " + str(round(Rloc * DIS_R[1], 2)) + " kpc" )
        result = plt.savefig('images/' + denFileName[8:-4] + '.png', dpi=200), print("go den/ to check " + denFileName[8:-4] + ".png")

    return result


def getFourierParasOf(denFileName, paraKey = 'a0'):
    '''
    HELPER function: get details about the fitting parameters
    Get the same parameter in an adjusted density file at each radii
    both of the two input are Strings.

    Arguments:
    denFileName: file stores number density to be Fourier fitted
    paraKey:     a0/a1/a2/a3/a4/a5/a6/w1/w2/w3/w4/w5/w6


    Returns:
    desired fitting parameter for denFileName
    '''
    valueBox = []

    for i in range(RESO_R):
        xdata = np.linspace(0, 2 * m.pi, RESO_THETA * 4)
        ydata = adjustFile(denFileName, 'y')[i]
        popt, pcov = curve_fit(FourierFunc, xdata, ydata)
        paraValues = np.abs(popt)
        paraKeys = ['a0','a1','a2','a3','a4','a5','a6','w1','w2','w3','w4','w5','w6']

        paraMap = dict(zip(paraKeys, paraValues))
        valueBox += [paraMap[paraKey]]

    return valueBox


def pitchFit(x, pConst, pK):
    result = pConst * np.exp(x * pK)
    return result

def getM2LogFit(denFileName, LEVEL=1.3):

    '''
    get the log fit for spirals with a mode of 2
    apologize for such a long function.
    (might be fixed in a non-near future :P

    Arguments:
    denFileName: file stores number density to be Fourier fitted
    LEVEL:       a threshold


    Returns:
    desired fitting parameter for denFileName
    '''

    OUTEDGE = int(14.5 / RM * RESO_R)
    INEDGE = int(4 / RM * RESO_R)

    CHOICE = input("Plot (Y/y) or Parameters (N/n): ") or "N"

    a2list = getFourierParasOf(denFileName, 'a2')
    a4list = getFourierParasOf(denFileName, 'a4')
    a6list = getFourierParasOf(denFileName, 'a6')

    mixList = np.dstack((a2list, a4list, a6list)).flatten()
    mixture = mixList.reshape(int(len(mixList) / 3), 3)

    indexArr = np.where(np.logical_and(LEVEL * mixture[:, 1] < mixture[:, 0], LEVEL * mixture[:, 2] < mixture[:, 0]))[0]
    listIndexArr = indexArr.tolist()
    indexStr = listIndexArr + [0]

    w2list = getFourierParasOf(denFileName, 'w2')

    if len(indexArr) == 0:
        result = print("Please reduce your LEVEL")

    else:
        boxA = []
        boxB = []

        while True:
            ruler = len(indexStr)
            if ruler <= 1:
                break
            boxB += [indexStr[0]]
            if indexStr[0] + 1 != indexStr[1]:
                boxA += [boxB]
                boxB = []
            indexStr.pop(0)

        lenOfSubList = [len(sublist) for sublist in boxA]
        wantedRange = boxA[lenOfSubList.index(max(lenOfSubList))]
        low = wantedRange[0]
        high = wantedRange[-1]

        if low < INEDGE:
            low = INEDGE
        if high > OUTEDGE:
            high = OUTEDGE


        xis = DIS_R[low: high]
        yis = w2list[low: high]

        i = 0
        while i < len(yis) - 1:
            if yis[i] < yis[i + 1] and (yis[i] + m.pi) > yis[i + 1]:
                i += 1
            else:
                if yis[i] > yis[i + 1]:
                    yis[i + 1] += m.pi
                else:
                    yis[i + 1] -= m.pi

            popt, pcov = curve_fit(pitchFit, xis, yis)
            pitchA = 90 - np.arctan(1 / popt[1]) * 180 / m.pi

            strPitch = str(round(pitchA, 3))

        if CHOICE == 'N' or CHOICE == 'n':
            print("-----FILE: " + denFileName + "-----")
            print("Estimated pitch angle( LAVEL = " + str(LEVEL) + ") is: " + strPitch + " degrees")
            print("Pconst is: " + str(round(popt[0], 4)))
            print("Pk is: " + str(round(popt[1], 4)))
            print("pcov is: " + str(pcov))
            result = print("----------Finished----------")

        else:
            (figure, axes) = plt.subplots()
            axes.yaxis.set_tick_params(direction='in', which='both')
            axes.xaxis.set_tick_params(direction='in', which='both')

            plt.scatter(xis, yis, label='data', color = COLOR2, alpha = 0.75)
            plt.plot(xis, pitchFit(xis, *popt), label='fit', linewidth=4, color = COLOR1)

            plt.xlabel('R [kpc]')
            plt.ylabel(r'$\phi$ [rad]')
            plt.title("Pitch Angle Log Fit \n " +  r'$\theta$ = ' + strPitch + r'$^{\circ}$')
            plt.legend(loc="lower right")
            result = plt.savefig('images/Fourier_' + denFileName[8:-4] + '.png', dpi=200), print("go den/ to check Fourier_" + denFileName[8:-4] + ".png")

    return result
