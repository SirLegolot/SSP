from    __future__ import division
from math import *
from visual import *

a = 2.773017979589484
e = 0.1750074901308245
M = radians(336.0050001501443)
Oprime = radians(108.032597191534)
iprime = radians(16.34548466739393)
wprime = radians(74.95130563682554)
sqrtmu = 0.01720209895
mu = sqrtmu**2
period = sqrt(4*pi**2*a**3/mu)
time = 0
deltat = 1

def solvekep(M):
    E = M # sets the intial guess of E to be equal to M
    Mguess = E - e*sin(E)
    while abs(Mguess - M) > 10**-4: # This for loop will loop until it finds a value for E that makes E-e*sinE very close to M
        E = E - (Mguess - M) / (1 - e*cos(E))
        Mguess = E - e*sin(E)
    return E

vecEcliptic = vector(0, 0, 0) # creates a vector for the asteroid, which will be the position of the asteroid
asteroid = sphere(pos=vecEcliptic*150, radius=(15), color=color.white) # this creates the asteroid
asteroid.trail = curve(color=color.white) # this sets up the trail of the asteroid
asteroid.label=label(pos=asteroid.pos, xoffset=20, yoffset=10, border=2, height=15, text='Asteroid') # this puts a label on asteroid
sun = sphere(pos=(0,0,0), radius=(50), color=color.yellow) # this creates the sun
sun.label=label(pos=sun.pos, xoffset=20, yoffset=10, border=2, height=15, text='Sun') # this puts a label on the sun

while True:
    rate(200)
    time = time + deltat
    Mcurrent = 2*pi/period*(time) + M # calculates the current M
    Ecurrent = solvekep(Mcurrent) # calculates the current E from M using the solve kep function
#Calculate x, y, and z coordinates
    cartesianX = a*(cos(Ecurrent)-e)
    cartesianY = a*sqrt(1-e**22)*sin(Ecurrent)
    cartesianZ = 0
#Calculates the radius of the asteroid's orbit
    radius = sqrt(cartesianX**2+cartesianY**2)
#Calculates the mean anomaly f
    cosf = cartesianX/radius
    sinf = cartesianY/radius
    f = atan2(sinf, cosf)
#Calculates the ecliptic coordinates of the asteroid
    eclipticX = radius*(cos(f+wprime)*cos(Oprime) - cos(iprime)*sin(f+wprime)*sin(Oprime))
    eclipticY = radius*(cos(iprime)*cos(Oprime)*sin(f+wprime) + cos(f+wprime)*sin(Oprime))
    eclipticZ = radius*sin(iprime)*sin(f+wprime)
#Sets the ecliptic coordinates as the x, y, and z values of the asteroid vector
    vecEcliptic = vector(eclipticX, eclipticY, eclipticZ)
#Updates the pos, trail, and label position
    asteroid.pos = vecEcliptic*150
    asteroid.trail.append(pos=asteroid.pos)
    asteroid.label.pos = asteroid.pos
