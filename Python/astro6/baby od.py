from __future__ import division
from math import *
import numpy as np

##vecrx = input("What is the x value of r vector: ")
##vecry = input("What is the y value of r vector: ")
##vecrz = input("What is the z value of r vector: ")
##vecrxdot = input("What is the x value of r dot vector: ")
##vecrydot = input("What is the y value of r dot vector: ")
##vecrzdot = input("What is the z value of r dot vector: ")

##vecr = np.array([vecrx, vecry, vecrz])
##vecrdot = np.array([vecrxdot, vecrydot, vecrydot])

def findQuadrant(sine, cosine):
    if cosine > 0 and sine > 0: #1
        return asin(sine)

    if cosine < 0 and sine > 0: #2
        return acos(cosine)

    if cosine < 0 and sine < 0: #3
        return pi - asin(sine)

    if cosine > 0 and sine < 0: #4
        return 2*pi + asin(sine)

mu = 0.01720209895**2
epsilon = radians(23.439281)

##vecr = np.array([0.244,2.17,-0.455])
##vecrdot = np.array([-0.731, -0.0041, 0.0502])
vecr = np.array([-1.481791,-1.52259774,0.56523413])
vecrdot = np.array([0.45254004, -0.12446764, -0.40393197])                

def babyOD(vecr, vecrdot):
    magvecr = sqrt(vecr[0]**2+vecr[1]**2+vecr[2]**2)


    a = ((2/magvecr)-np.dot(vecrdot, vecrdot))**-1
    print "a: ",a

    vech = np.cross(vecr, vecrdot)
    magvech = sqrt(vech[0]**2+vech[1]**2+vech[2]**2)
    e = sqrt(1-magvech**2/a)
    print "e: ", e

    i = acos(vech[2]/magvech)
    print "i: ", degrees(i)

    sinO = vech[0]/(magvech*sin(i))
    cosO = -vech[1]/(magvech*sin(i))
    O = findQuadrant(sinO, cosO)
    print "O: ", degrees(O)

    sinU = vecr[2]/(magvecr*sin(i))
    cosU = (vecr[0]*cos(O)+vecr[1]*sin(O))/magvecr
    U = findQuadrant(sinU, cosU)

    sinV = (a*(1-e**2)*np.dot(vecr, vecrdot))/(magvech*magvecr*e)
    cosV = ((a*(1-e**2)/magvecr)-1)/e
    V = findQuadrant(sinV, cosV)

    w = U-V
    if w<0:
        print "w: ", degrees(w)+360
    else:
        print "w: ", degrees(w)

    E = acos((1-magvecr/a)*(1/e))
    M = E - e*sin(E)
    if V>0:
        print "M: ", 360-degrees(M)
    else:
        print "M: ", M

    return a, e, i, O, w, M

print babyOD(vecr, vecrdot)





