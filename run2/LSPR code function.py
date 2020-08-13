from __future__ import division
from math import *
import numpy as np
from numpy.linalg import inv 

astrdCentroidx = input("What is the unknown object's x centroid value: ")
astrdCentroidy = input("What is the unknown object's y centroid value: ")

def lspr(filename, centroidx, centroidy):
    file = np.loadtxt(filename)
    ##print "file: "
    ##print file

    astrdCentroidx = centroidx
    astrdCentroidy = centroidy


    sumAi = np.sum(file[:,2])
    ##print "sumAi: ", sumAi

    listAiXi = []
    for i in file[:,::2]:
        listAiXi.append(i[0]*i[1])
    sumAiXi = np.sum(listAiXi)
    ##print "sumAiXi: ", sumAiXi

    listAiYi = []
    for i in file[:,1:3]:
        listAiYi.append(i[0]*i[1])
    sumAiYi = np.sum(listAiYi)
    ##print "sumAiYi: ", sumAiYi

    N = len(file)
    ##print "N: ", N

    sumXi = np.sum(file[:,0])
    ##print "sumXi: ", sumXi

    sumYi = np.sum(file[:,1])
    ##print "sumYi: ", sumYi

    listXi2 = []
    for i in file[:,0]:
        listXi2.append(i**2)
    sumXi2 = np.sum(listXi2)
    ##print "sumXi2: ", sumXi2

    listXiYi = []
    for i in file[:,:2]:
        listXiYi.append(i[0]*i[1])
    sumXiYi = np.sum(listXiYi)
    ##print "sumXiYi: ", sumXiYi

    listYi2 = []
    for i in file[:,1]:
        listYi2.append(i**2)
    sumYi2 = np.sum(listYi2)
    ##print "sumYi2: ", sumYi2

    raMatrix = np.array([[sumAi],
                        [sumAiXi],
                        [sumAiYi]])
    ##print "raMatrix: "
    ##print raMatrix

    raMatrixToBeInv = np.array([[N, sumXi, sumYi],
                                [sumXi, sumXi2, sumXiYi],
                                [sumYi, sumXiYi, sumYi2]])
    ##print "raMatrixToBeInv: "
    ##print raMatrixToBeInv

    raMatrixInv = np.linalg.inv(raMatrixToBeInv)
    ##print "raMatrixInv: "
    ##print raMatrixInv

    b1a11a12 = np.array([[sumAi*raMatrixInv[0,0]+sumAiXi*raMatrixInv[0,1]+sumAiYi*raMatrixInv[0,2]],
                         [sumAi*raMatrixInv[1,0]+sumAiXi*raMatrixInv[1,1]+sumAiYi*raMatrixInv[1,2]],
                         [sumAi*raMatrixInv[2,0]+sumAiXi*raMatrixInv[2,1]+sumAiYi*raMatrixInv[2,2]]])
    ##print "b1a11a12:"
    ##print b1a11a12

    b1 = b1a11a12[0,0]
    ##print "b1: ", b1

    a11 = b1a11a12[1,0]
    ##print "a11: ", a11

    a12 = b1a11a12[2,0]
    ##print "a12: ", a12

    sumDi = np.sum(file[:,3])
    ##print "sumDi: ", sumDi

    listDiXi = []
    for i in file[:,::3]:
        listDiXi.append(i[0]*i[1])
    sumDiXi = np.sum(listDiXi)
    ##print "sumDiXi: ", sumDiXi

    listDiYi = []
    for i in file[:,1::2]:
        listDiYi.append(i[0]*i[1])
    sumDiYi = np.sum(listDiYi)
    ##print "sumDiYi: ", sumDiYi

    b2a21a22 = np.array([[sumDi*raMatrixInv[0,0]+sumDiXi*raMatrixInv[0,1]+sumDiYi*raMatrixInv[0,2]],
                         [sumDi*raMatrixInv[1,0]+sumDiXi*raMatrixInv[1,1]+sumDiYi*raMatrixInv[1,2]],
                         [sumDi*raMatrixInv[2,0]+sumDiXi*raMatrixInv[2,1]+sumDiYi*raMatrixInv[2,2]]])
    ##print "b2a21a22:"
    ##print b2a21a22

    b2 = b2a21a22[0,0]
    ##print "b2: ", b2

    a21 = b2a21a22[1,0]
    ##print "a21: ", a21

    a22 = b2a21a22[2,0]
    ##print "a22: ", a22

    listsumuncrtnA = []
    for i in range(N):
        listsumuncrtnA.append((file[i, 2]-b1-a11*file[i,0]-a12*file[i,1])**2)
    sumuncrtnA = np.sum(listsumuncrtnA)
    uncrtnA = sqrt(sumuncrtnA/(N-3))*3600
    printUncrtnA = '%.2f' % round(uncrtnA,2)
    ##print "uncrtnA: ", uncrtnA

    listsumuncrtnD = []
    for i in range(N):
        listsumuncrtnD.append((file[i, 3]-b2-a21*file[i,0]-a22*file[i,1])**2)
    sumuncrtnD = np.sum(listsumuncrtnD)
    uncrtnD = sqrt(sumuncrtnD/(N-3))*3600
    printUncrtnD = '%.2f' % round(uncrtnD,2)
    ##print "uncrtnD: ", uncrtnD

    RA = (b1+a11*astrdCentroidx+a12*astrdCentroidy)/15
    RAh = int(RA)
    RAmin = int((RA-RAh)*60)
    RAsec = (RA-RAh-(RAmin/60))*3600
    printRAsec = '%.2f' % round(RAsec,2)

    Dec = b2+a21*astrdCentroidx+a22*astrdCentroidy
    absDec = abs(Dec)
    Decdeg = int(absDec)
    Decarcmin = int((absDec-Decdeg)*60)
    Decarcsec = ((absDec-Decdeg)-(Decarcmin/60))*3600
    printDecarcsec = '%.2f' % round(Decarcsec,2)

    return b1, b2, a11, a12, a21, a22, printUncrtnA, printUncrtnD, astrdCentroidx, astrdCentroidy, RAh, RAmin, printRAsec, Dec, Decdeg, Decarcmin, printDecarcsec

b1, b2, a11, a12, a21, a22, printUncrtnA, printUncrtnD, astrdCentroidx, astrdCentroidy, RAh, RAmin, printRAsec, Dec, Decdeg, Decarcmin, printDecarcsec = lspr("astrom2a.txt", astrdCentroidx,  astrdCentroidy)

print "****************"
print "plate constants:"
print "****************"
print "b1: ", b1, "deg"
print "b2: ", b2, "deg"
print "a11: ", a11, "deg/pix"
print "a12: ", a12, "deg/pix"
print "a21: ", a21, "deg/pix"
print "a22: ", a22, "deg/pix"
print "***********"
print "uncertainty"
print "***********"
print "uncrtnA: ", printUncrtnA, "arcsec"
print "uncrtnD: ", printUncrtnD, "arcsec"
print "***********************************"
print "astrometry for (x,y)=( "+ str(astrdCentroidx)+ ", "+ str(astrdCentroidy)+ ")"
print "***********************************"
print "RA =", str(RAh)+":"+str(RAmin)+":"+printRAsec
if Dec>=0:
    print "Dec =", "+"+str(Decdeg)+":"+str(Decarcmin)+":"+printDecarcsec
else:
    print "Dec =", "-"+str(Decdeg)+":"+str(Decarcmin)+":"+printDecarcsec
