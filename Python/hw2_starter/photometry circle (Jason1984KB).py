from __future__ import division
import numpy as np
import matplotlib.pyplot as plt
import math
from astropy.io import fits

myim = fits.getdata("Average1.fit")
##plt.imshow(myim, vmin=myim.mean(), vmax=2*myim.mean())
##plt.gray()
##plt.show()


##starx = input("What is your star x coordinate?" )
##stary = input("What is your star y coordinate?" )
##magstar = input("What is your star's magnitude" )
##blankx = input("What is your blank x coordinate?" )
##blanky = input("What is your blank y coordinate?" )
##astrdx = input("What is your asteroid x coordinate?" )
##astrdy = input("What is your asteroid y coordinate?" )

#Jason's data:

##myim
##radius =7,6

starx = 866
stary = 164
magstar = 13.66
blankx = 63
blanky = 910
astrdx = 584
astrdy = 549


print "star box"
starbox = myim[stary-6:stary+7, starx-6:starx+7]#extracts a 14 by 14 grid around the given star coordinates
y,x =np.ogrid[-6:7,-6:7]
mask = x**2 +y**2 >= 7**2
starbox[mask] = 0
print starbox

print "total star (star+sky) " # gets the total of all the pixels in the star box
totalstar = sum(sum(starbox))
print totalstar








print "blank box"
blankbox = myim[blanky-6:blanky+7, blankx-6:blankx+7] #extracts a 14 by 14 grid around the given coordinates
blankbox[mask] = 0
print blankbox

print "total blank " # gets the total of all the pixels in the blank box
totalblank = sum(sum(blankbox))
print totalblank

print "average blank (avgSky)" # gets the average pixel value in the blank box
averageblank = totalblank/((13**2)-24)
print averageblank

print "signal star" # calculates the signal for star
signal = totalstar - (averageblank*((13**2)-24))
print signal

print "constant" # gets the value of constant
const = magstar + 2.5*math.log10(signal)
print const







print "asteroid box"
astrdbox = myim[astrdy-6:astrdy+7, astrdx-6:astrdx+7] #extracts a 14 by 14 grid around the given asteroid coordinates
astrdbox[mask] = 0
print astrdbox

print "total astrd (star+sky) " # gets the total of all the pixels in the asteroid box
totalastrd = sum(sum(astrdbox))
print totalastrd


print "signal asteroid" # calculates the signal of asteroid
signalastrd = totalastrd - (averageblank*((13**2)-24))
print signalastrd

print "magnitude asteroid" # calculates magnitude of the asteroid
magastrd = -2.5*math.log10(signalastrd) + const
print magastrd
