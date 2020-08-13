from    __future__ import division
from math import *
import numpy as np



def ephem(a, e, i, O, w, M0, t, t0, earthSun):
    epsilon = radians(23.439281)
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
    radius = sqrt(cartesian[0]**2+ cartesian[1]**2)

    #Calculates the mean anomaly f
    cosf = cartesian[0]/radius
    sinf = cartesian[1]/radius
    f = atan2(sinf, cosf)
    if f<0:
        f=f+2*pi

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

    #RA
    sinRA = equatorial[1]/(magRangeVec*cos(Dec))
    cosRA = equatorial[0]/(magRangeVec*cos(Dec))
    RA = atan2(sinRA, cosRA)
    if RA<0:
        RA = RA+2*pi

    return RA, Dec




