import seaborn as sns
import os
import math as m
import matplotlib.pyplot as plt
from matplotlib.pyplot import figure
from numpy import loadtxt
from scipy.optimize import curve_fit
from scipy import optimize
import numpy as np
from phycon import *
from grid import *

#-------------------------------- Fourier --------------------------------------#

# model 
def FourierFunc(x, a0, a1, a2, a3, a4, a5, a6, w1, w2, w3, w4, w5, w6):
    '''
    Fourier serise model for fitting 
    '''
    part1 = a0 + a1 * np.cos(x - w1) + a2 * np.cos(2 * (x - w2)) + a3 * np.cos(3 * (x - w3))
    part2 = a4 * np.cos(4 * (x - w4)) + a5 * np.cos(5 * (x - w5)) + a6 * np.cos(6 * (x - w6))
    result = part1 + part2
    return result

# Fitting process --- for plotting
def fourierFit(denFileName, Rloc = 40, parameters = 'y'):
    '''
    Fit density series at each radius with:
    'y' for printing parameters dictionary
    'n' for plotting fitting curve
    
    NOTE: adjusted_denfile is one line in that file instead of the whole
    '''
    xdata = np.linspace(0, 2 * m.pi, RESO_THETA * 4) 
    ydata = adjustFile(denFileName)[Rloc]

    popt, pcov = curve_fit(FourierFunc, xdata, ydata)
    
    if parameters == 'y' or parameters == 'Y':
        paraValues = np.abs(popt)
        paraKeys = ['a0','a1','a2','a3','a4','a5','a6','w1','w2','w3','w4','w5','w6']
        result = dict(zip(paraKeys, paraValues))
        # print("----------------------------")
        # print("Fourier fit at R = " + str(round(Rloc * DIS_R[1], 2)) + " kpc is :")
        # print(result)
        # print("----------------------------")
    
    if parameters == 'n' or parameters == 'N':
        
        plt.plot(xdata, ydata, label='data')
        plt.plot(xdata, FourierFunc(xdata, *popt),  label='fit', linewidth = 4)
        plt.xlabel('x [kpc]')
        plt.ylabel('y [kpc]')
        plt.legend(loc = "lower right")
        plt.title("Fourier fit at R = " + str(round(Rloc * DIS_R[1], 2)) + " kpc")
        result = plt.show()
    
    return result

# get details about the fitting parameters
def getFourierParasOf(denFileName, paraKey):
    '''
    Get the same paramter in an adjusted density file at each redii
    both of the two input are in the form of string.
    '''
    valueBox = []
    for i in range(RESO_R):
        paraMap = fourierFit(denFileName, Rloc = i, parameters = 'y')
        valueBox += [paraMap[paraKey]]
    
    #print(np.array(valueBox).reshape(4, 20))
        
    return valueBox


def pitchFit(x, pConst, pK):
    result = pConst * np.exp(x * pK)
    return result

def getM2LogFit(denFileName, plot = 'Y', LEVEL = 1.3):
    w2list = getFourierParasOf(denFileName, 'w2')
    w4list = getFourierParasOf(denFileName, 'w4')
    w6list = getFourierParasOf(denFileName, 'w6')
    
    mixList = np.dstack((w2list, w4list, w6list)).flatten()
    mixture = mixList.reshape(int(len(mixList) / 3), 3)
    
    indexArr = np.where(np.logical_and(LEVEL * mixture[:, 1] < mixture[:, 0], LEVEL * mixture[:, 2] < mixture[:, 0]))[0]
    indexStr = indexArr.tolist() + [0]
    
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
        high = wantedRange[-1] - 4 #! !boundary effect 
    
        xis = DIS_R[low : high]
    
        yis = w2list[low : high]
        i = 0
        while i < len(yis) - 1:
            if yis[i] < yis[i + 1] and (yis[i] + m.pi) > yis[i + 1] :
                i += 1
            else:
                if yis[i] > yis[i + 1]:
                    yis[i + 1] += m.pi
                else:
                    yis[i + 1] -= m.pi            

            popt, pcov = curve_fit(pitchFit, xis, yis) 
            pitchA = 90 - np.arctan(1 / popt[1]) * 180 / m.pi
            strPitch = str(round(pitchA,3))
        
        if plot == 'N' or plot == 'n':
            print("-----FILE: " + denFileName + "-----")
            print("Estimated pitch angle( LAVEL = " + str(LEVEL) + ") is: " + strPitch + " degrees")
            print("Pconst is: " + str(round(popt[0],4)))
            print("Pk is: " + str(round(popt[1],4)))
            print("pcov is: " + str(pcov))
            result = print("----------Finished----------")
        
        else:
            plt.scatter(xis, yis, label='data')
            plt.plot(xis, pitchFit(xis, *popt),  label='fit', linewidth = 4)

            plt.xlabel('R [kpc]')
            plt.ylabel(r'$\phi$')
            plt.legend(loc = "lower right")
            plt.title("Logarithmic fit for mode = 2")
            #plt.text(11.5, 35, r"$\theta_{p}=$" + strPitch, bbox={'facecolor': 'white', 'pad': 2.7})
            result = plt.show() 
    
    return result 