from __future__ import division
import numpy as np
import matplotlib.pyplot as plt
import math
from astropy.io import fits


im = fits.getdata("sampleimage.fits")
myim = fits.getdata("Average1.fit")
#plt.imshow(im, vmin=im.mean(), vmax=2*im.mean())
#plt.gray()
#plt.show()


##starx = input("What is your star x coordinate?" )
##stary = input("What is your star y coordinate?" )
##magstar = input("What is your star's magnitude" )
##blankx = input("What is your blank x coordinate?" )
##blanky = input("What is your blank y coordinate?" )
##astrdx = input("What is your asteroid x coordinate?" )
##astrdy = input("What is your asteroid y coordinate?" )

#Example data
starx = 173
stary = 342
magstar = 15.26
blankx = 219
blanky = 163
astrdx = 351
astrdy = 154

#Jason's data:

##myim
##radius =7,6
##starx = 866
##stary = 164
##magstar = 13.66
##blankx = 63
##blanky = 910
##astrdx = 584
##astrdy = 549


print "star box"
starbox = im[stary-5:stary+4, starx-4:starx+5]#extracts a 14 by 14 grid around the given star coordinates
print starbox
y,x =np.ogrid[-4:5,-4:5]
mask = x**2 +y**2 >= 5**2
starbox[mask] = 0
print starbox

print "total star (star+sky) " # gets the total of all the pixels in the star box
totalstar = sum(sum(starbox))
print totalstar








print "blank box"
blankbox = im[blanky-5:blanky+4, blankx-4:blankx+5] #extracts a 14 by 14 grid around the given coordinates
blankbox[mask] = 0
print blankbox

print "total blank " # gets the total of all the pixels in the blank box
totalblank = sum(sum(blankbox))
print totalblank

print "average blank (avgSky)" # gets the average pixel value in the blank box
averageblank = totalblank/(9**2)
print averageblank

print "signal star" # calculates the signal for star
signal = totalstar - (averageblank*(9**2))
print signal

print "constant" # gets the value of constant
const = magstar + 2.5*math.log10(signal)
print const







print "asteroid box"
astrdbox = im[astrdy-5:astrdy+4, astrdx-4:astrdx+5] #extracts a 14 by 14 grid around the given asteroid coordinates
astrdbox[mask] = 0
print astrdbox

print "total astrd (star+sky) " # gets the total of all the pixels in the asteroid box
totalastrd = sum(sum(astrdbox))
print totalastrd


print "signal asteroid" # calculates the signal of asteroid
signalastrd = totalastrd - (averageblank*(9**2))
print signalastrd

print "magnitude asteroid" # calculates magnitude of the asteroid
magastrd = -2.5*math.log10(signalastrd) + const
print magastrd