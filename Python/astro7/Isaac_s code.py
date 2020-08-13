from __future__ import division
import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
from math import *

input_image = fits.getdata('run 1.fit')
starAptFull = np.copy(input_image)
astroAptFull = np.copy(input_image)
blankAnnulusFull = np.copy(input_image)

coorx= 767 # input("Enter the Star X coordinates:") 
coory= 314 # input("Enter the Star Y coordinates:")
mag = 14.54
astroidx= 442 # input("Enter the Asteroid X coordinates:") 
astroidy= 473 # input("Enter the Asteroid Y coordinates:")
rCntrd= 2 #input("Enter the width (preferably an odd number:")

rAstro= 6

rApt= 6 #input("Enter the desired radius of the circular aperture:")
rIn= 9
rEx= 12
buf =5

median= np.median(input_image)

starcentroid= input_image-median
astrocentroid= input_image-median

matrix = starcentroid[coory-rCntrd:coory+rCntrd+1,coorx-rCntrd:coorx+rCntrd+1]
print matrix
#matrix = input_image[coory-width:coory+width+1,coorx-width:coorx+width+1]

matrix2 = astrocentroid[astroidy-rCntrd:astroidy+rCntrd+1,astroidx-rCntrd:astroidx+rCntrd+1]
print matrix2

alltotal= matrix.sum() # The total value after summing up the entire matrix
alltotal2= matrix2.sum()

xtotal= sum(matrix)
xtotal2= sum(matrix2)

ytotal=[]
for i in matrix:    #Gives out a list of the sums of each row of data
    ytotal.append(sum(i))
    
ytotal2=[]
for i in matrix2:    #Gives out a list of the sums of each row of data
    ytotal2.append(sum(i))

def centroidcalculationstar(value,coor): # Calculation 
    centroid = 0
    for counter in range(len(value)):
        centroid= centroid + (value[counter]*counter)
    finalcentroid = (centroid/(alltotal))
    finalfinal= finalcentroid+ coor -1

    return finalfinal

def centroidcalculationastro(value,coor): # Calculation 
    centroid = 0
    for counter in range(len(value)):
        centroid= centroid + (value[counter]*counter)
    finalcentroid = (centroid/(alltotal2))
    finalfinal= finalcentroid+ coor -1

    return finalfinal

xstar= centroidcalculationstar(xtotal, coorx)
ystar= centroidcalculationstar(ytotal, coory)
xastro= centroidcalculationastro(xtotal2, astroidx)
yastro= centroidcalculationastro(ytotal2, astroidy)

print xstar, ystar, xastro, yastro

deltaStar_X= coorx - xstar
deltaStar_Y= coory - ystar

deltaAstro_X= astroidx - xastro
deltaAstro_Y= astroidy - yastro

#############################################
starApt = starAptFull[coory-rApt-buf:coory +rApt+buf+1,coorx-rApt-buf:coorx+buf+rApt+1]

y,x = np.ogrid[-rApt-buf:rApt+buf+1, -rApt-buf:rApt+buf+1] #circular aperture
mask1 = (x+deltaStar_X)**2 + (y+deltaStar_Y)**2 >= rApt**2

starApt[mask1]=0
print starApt

#############################################
blankAnnulus = blankAnnulusFull[coory-rEx-buf:coory +buf+rEx+1,coorx-buf-rEx:coorx+buf+rEx+1]

y2,x2 = np.ogrid[-rEx-buf:rEx+buf+1, -rEx-buf:rEx+buf+1] #circular aperture
mask2 = (x2+deltaStar_X)**2 + (y2+deltaStar_Y)**2 >= rEx**2
mask3 = (x2+deltaStar_X)**2 + (y2+deltaStar_Y)**2 <= rIn**2

blankAnnulus[mask2]=0
blankAnnulus[mask3]=0

print blankAnnulus
#######################   Asteroid   ########################### 
astroApt = astroAptFull[astroidy-rApt-buf:astroidy +buf+rApt+1,astroidx-buf-rApt:astroidx+buf+rApt+1]

y3,x3 = np.ogrid[-rAstro-buf:rAstro+buf+1, -rAstro-buf:rAstro+buf+1] #circular aperture # To change, change rApt to rAstro
mask4 = (x+deltaStar_X)**2 + (y+deltaStar_Y)**2 >= rAstro**2 # Can change if there's a need for a difference in radii

astroApt[mask4]=0
print astroApt

totalADUAnnulus= np.sum(blankAnnulus)
totalADUStar= np.sum(starApt)

filledPixels_Annulus= np.count_nonzero(blankAnnulus)
filledPixels_Star= np.count_nonzero(starApt)

avgSky= totalADUAnnulus/filledPixels_Annulus
print "avgSky", avgSky
signal= totalADUStar-(avgSky*(filledPixels_Star))
print "signal",signal

constant= mag+(2.5*(log10(signal)))
print "constant", constant

totalADUastro= np.sum(astroApt)
filledPixels_Astro= np.count_nonzero(astroApt)

astroSignal= totalADUastro-(avgSky*(filledPixels_Astro))

print "astroSignal",astroSignal

astromag= (-2.5*(log10(astroSignal)))+constant
print astromag
