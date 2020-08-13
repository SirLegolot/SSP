from __future__ import division
import numpy as np
from visual import *
from math import *

def ephem(a,e,i,O,w,Mo,t,to,EarthSunVec):
    sqrtmu1 = 0.01720209895
    mu1 = sqrtmu1**2
    epsilon= radians(23.439281)

    n = sqrt(mu1)/sqrt(a**3)

    M = Mo + (n*(t-to))

    def solvekep(M):
        Eguess = M
        Mguess = Eguess - e*sin(Eguess)
        while abs(Mguess - M) > 1e-004:
            Mguess = Eguess - e*sin(Eguess)
            Eguess = Eguess - (Eguess - e*sin(Eguess) - M) / (1 - e*cos(Eguess))
        return Eguess

    Etrue = solvekep(M)

    cartesianx= a*((cos(Etrue)-e))
    cartesiany= a*sin(Etrue)*(sqrt(1-(e**2)))
    cartesianz= 0

    radius= sqrt((cartesianx**2)+(cartesiany**2))
    cosinef= cartesianx/radius
    sinef= cartesiany/radius
    f= atan2(sinef,cosinef)

    ecliptic_x = radius*(((cos(f+w)*cos(O)) - (sin(f+w)*cos(i)*sin(O))))
    ecliptic_y = radius*(((cos(i)*cos(O)*sin(f+w))) + ((cos(f+w)*sin(O))))
    ecliptic_z = radius*sin(i)*sin(f+w)

    #### Earth to Sun Position Vectors ########

    EarthSunXVec= EarthSunVec[0]
    EarthSunYVec= EarthSunVec[1]
    EarthSunZVec= EarthSunVec[2]

    ####### Earth to Asteroid ##########

    EarthAsteroidXVec =  ecliptic_x+EarthSunXVec
    EarthAsteroidYVec = ecliptic_y+EarthSunYVec
    EarthAsteroidZVec = ecliptic_z+EarthSunZVec

    mag = sqrt((EarthAsteroidXVec**2)+(EarthAsteroidYVec**2)+(EarthAsteroidZVec**2))

    ####### Earth to Asteroid in Equatorial Coordinates #########

    EAEquatorialX= EarthAsteroidXVec
    EAEquatorialY= (EarthAsteroidYVec*cos(epsilon))-(EarthAsteroidZVec*sin(epsilon))
    EAEquatorialZ= (EarthAsteroidYVec*sin(epsilon))+(EarthAsteroidZVec*cos(epsilon))

    ####### Converting to RA and DEC ##############

    delta= asin(EAEquatorialZ/mag)
    cosinealpha= EAEquatorialX/(mag*cos(radians(delta)))     # cosine of an angle in radians goes into atan2 
    sinealpha= EAEquatorialY/(mag*cos(radians(delta)))       #this is sine of the same angle in radians
    alpha= atan2(sinealpha,cosinealpha)             #actual angle

        
    

    return alpha, delta
