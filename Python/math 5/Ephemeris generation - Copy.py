from    __future__ import division
from math import *
import numpy as np

e = 0.4518788622141914
a = 1.704987880150487
i = radians(4.620927383463039)
O = radians(144.567249166763)
w = radians(197.6954745176973)
M0 = radians(287.3388829192644)
t = 2456842.5
t0 = 2456000.5
epsilon = radians(23.4347)
earthSun = np.array([-2.074064555488755*10**-1, 9.953010074087935*10**-1, -3.313757729661123*10**-5])
##e = 0.09618818215877029
##a = 1.844195082375194
##i = radians(23.66084024789276)
##O = radians(132.0979894264106)
##w = radians(129.1630964507269)
##M0 = radians(113.1648904822557)

def ephem(a, e, i, O, w, M0, t0, t, earthSun):
    mu = 0.01720209895**2
    
    #calculate n
    n = sqrt(mu/a**3)

    #calculate mean anomaly
    M = M0+n*(t-t0)
    
    #Newton's method to calculate eccentric anomaly
    def solvekep(Mtrue):
        E = Mtrue # sets the intial guess of E to be equal to M
        Mguess = E - e*sin(E)
        while abs(Mguess - Mtrue) > 10**-4: # This for loop will loop until it finds a value for E that makes E-e*sinE very close to M
            Mguess = E - e*sin(E)
            E = E - (Mguess - Mtrue) / (1 - e*cos(E))
        return E
    E = solvekep(M)
    
    #Calculate x, y, and z cartesian coordinates
    cartesian = np.array([a*(cos(E)-e),a*sqrt(1-e**2)*sin(E), 0])


    #Calculates the radius of the asteroid's orbit
    radius = np.linalg.norm(cartesian)

    #Calculates the mean anomaly f
    cosf = cartesian[0]/radius
    sinf = cartesian[1]/radius
    f = atan2(sinf, cosf)

    #Calculates the ecliptic coordinates of the asteroid
    ecliptic = np.array([radius*(cos(f+w)*cos(O) - cos(i)*sin(f+w)*sin(O)), radius*(cos(i)*cos(O)*sin(f+w) + cos(f+w)*sin(O)), radius*sin(i)*sin(f+w)])


    #range vector
    rangeVec = np.array([ecliptic[0] + earthSun[0], ecliptic[1] + earthSun[1], ecliptic[2] + earthSun[2]])
    magRangeVec = np.linalg.norm(rangeVec)


    #range unit vector
    rangeUnitVec = rangeVec/magRangeVec

    # equatorial coordinates
    equatorial = np.array([rangeVec[0], rangeVec[1]*cos(epsilon) - rangeVec[2]*sin(epsilon), rangeVec[1]*sin(epsilon) + rangeVec[2]*cos(epsilon)])

    #Dec
    Dec = asin(equatorial[2]/magRangeVec)
    degDec = degrees(Dec)
    absDec = abs(degDec)
    Decdeg = int(absDec)
    Decarcmin = int((absDec-Decdeg)*60)
    Decarcsec = ((absDec-Decdeg)-(Decarcmin/60))*3600
    roundDecarcsec = float('%.1f' % round(Decarcsec,2))
    if Dec>=0:
        Decfinal = radians(Decdeg + Decarcmin/60 + roundDecarcsec/3600)
    else:
        Decfinal = radians(-(Decdeg + Decarcmin/60 + roundDecarcsec/3600)) 

    #RA
    sinRA = equatorial[1]/(magRangeVec*cos(Dec))
    cosRA = equatorial[0]/(magRangeVec*cos(Dec))
    RA = atan2(sinRA, cosRA)
    if RA<0:
        degRA = (degrees(RA)+360)/15
    else:
        degRA = degrees(RA)/15
    RAh = int(degRA)
    RAmin = int((degRA-RAh)*60)
    RAsec = (degRA-RAh-(RAmin/60))*3600
    roundRAsec = float('%.2f' % round(RAsec,2))
    RAfinal = radians((RAh + RAmin/60 + roundRAsec/3600)*15)

    return RAfinal, Decfinal

print ephem(a, e, i, O, w, M0, t0, t, earthSun)



