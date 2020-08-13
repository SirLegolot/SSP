from    __future__ import division
from math import *
import numpy as np

e = 0.4518788622141914
a = 1.704987880150487
i = radians(4.620927383463039)
O = radians(144.567249166763)
w = radians(197.6954745176973)
M0 = radians(287.3388829192644)
##e = 0.09618818215877029
##a = 1.844195082375194
##i = radians(23.66084024789276)
##O = radians(132.0979894264106)
##w = radians(129.1630964507269)
##M0 = radians(113.1648904822557)

def ephem(a, e, i, O, w, M0, t0, t)
    mu = 0.01720209895**2
    t = 2456842.5
    t0 = 2456000.5
    ##t = 2458118.5
    ##t0 = 2458000.5
    epsilon = radians(23.4347)
    #calculate n
    n = sqrt(mu/a**3)

    #calculate mean anomaly
    M = M0+n*(t-t0)
    ##print M
    #Newton's method to calculate eccentric anomaly
    def solvekep(Mtrue):
        E = Mtrue # sets the intial guess of E to be equal to M
        Mguess = E - e*sin(E)
        while abs(Mguess - Mtrue) > 10**-4: # This for loop will loop until it finds a value for E that makes E-e*sinE very close to M
            Mguess = E - e*sin(E)
            E = E - (Mguess - Mtrue) / (1 - e*cos(E))
        return E
    E = solvekep(M)
    ##print E
    #Calculate x, y, and z cartesian coordinates
    cartesianX = a*(cos(E)-e)
    cartesianY = a*sqrt(1-e**2)*sin(E)
    cartesianZ = 0
    ##print cartesianX, cartesianY, cartesianZ


    #Calculates the radius of the asteroid's orbit
    radius = sqrt(cartesianX**2+cartesianY**2)

    #Calculates the mean anomaly f
    cosf = cartesianX/radius
    sinf = cartesianY/radius
    f = atan2(sinf, cosf)

    #Calculates the ecliptic coordinates of the asteroid
    eclipticX = radius*(cos(f+w)*cos(O) - cos(i)*sin(f+w)*sin(O))
    eclipticY = radius*(cos(i)*cos(O)*sin(f+w) + cos(f+w)*sin(O))
    eclipticZ = radius*sin(i)*sin(f+w)
    ##print eclipticX, eclipticY, eclipticZ

    #Earth-sun vector
    earthSunX = -2.074064555488755*10**-1
    earthSunY = 9.953010074087935*10**-1
    earthSunZ = -3.313757729661123*10**-5
    ##earthSunX = 1.579922623076506*10**-1
    ##earthSunY = -9.705438166221750*10**-1
    ##earthSunZ = 3.780898978569830*10**-5

    #range vector
    rangeVecX = eclipticX + earthSunX
    rangeVecY = eclipticY + earthSunY
    rangeVecZ = eclipticZ + earthSunZ
    magRangeVec = sqrt(rangeVecX**2+rangeVecY**2+rangeVecZ**2)
    ##print magRangeVec
    ##print rangeVecX, rangeVecY, rangeVecZ

    #range unit vector
    rangeUnitVecX = rangeVecX/magRangeVec
    rangeUnitVecY = rangeVecY/magRangeVec
    rangeUnitVecZ = rangeVecZ/magRangeVec
    ##print rangeUnitVecX, rangeUnitVecY, rangeUnitVecZ

    # equatorial coordinates
    equatorialX = rangeVecX
    equatorialY = rangeVecY*cos(epsilon) - rangeVecZ*sin(epsilon)
    equatorialZ = rangeVecY*sin(epsilon) + rangeVecZ*cos(epsilon)
    ##print equatorialX, equatorialY, equatorialZ

    #Dec
    Dec = asin(equatorialZ/magRangeVec)
    degDec = degrees(Dec)
    absDec = abs(degDec)
    Decdeg = int(absDec)
    Decarcmin = int((absDec-Decdeg)*60)
    Decarcsec = ((absDec-Decdeg)-(Decarcmin/60))*3600
    printDecarcsec = '%.2f' % round(Decarcsec,2)

#RA
sinRA = equatorialY/(magRangeVec*cos(Dec))
cosRA = equatorialX/(magRangeVec*cos(Dec))
RA = atan2(sinRA, cosRA)
if RA<0:
    degRA = (degrees(RA)+360)/15
else:
    degRA = degrees(RA)/15
RAh = int(degRA)
RAmin = int((degRA-RAh)*60)
RAsec = (degRA-RAh-(RAmin/60))*3600
printRAsec = '%.2f' % round(RAsec,2)

print "RA =", str(RAh)+":"+str(RAmin)+":"+printRAsec
if Dec>=0:
    print "Dec =", "+"+str(Decdeg)+":"+str(Decarcmin)+":"+printDecarcsec
else:
    print "Dec =", "-"+str(Decdeg)+":"+str(Decarcmin)+":"+printDecarcsec



