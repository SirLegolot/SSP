from __future__ import division
import numpy as np
import matplotlib.pyplot as plt
import math
from astropy.io import fits

myim = fits.getdata("run 1.fit")
myimcentrd = np.copy(myim)-np.median(myim)
myimstar = np.copy(myim)
myimblank = np.copy(myim)
myimastrd = np.copy(myim)
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

starx = 767
stary = 314
magstar = 14.54
astrdx = 442
astrdy = 473
radiusc = 2
radius = 6
radius2 = 9
radius3 = 12
buf = 3

def calculateCentroid(xcoord, ycoord, radius):
    #Creates nine point array centered on the coordinates
    data = myimcentrd[ycoord-radius:ycoord+radius+1, xcoord-radius:xcoord+radius+1]

    # totals list x
    totalListx = 0
    for i in data[:,:]:
        totalListx = totalListx + i

    # weighted total x
    wghtTotalx = 0
    for index in range(len(totalListx)):
        wghtTotalx = totalListx[index]*index + wghtTotalx

    # centroid x
    centroidx = wghtTotalx/data.sum()
    centroidxdisplay = centroidx +xcoord-radius


    # totals list y
    totalListy = []
    for row in data:
        totalListy.append(np.sum(row))

    # weighted total y
    wghtTotaly = 0
    for index in range(len(totalListy)):
        wghtTotaly = totalListy[index]*index + wghtTotaly

    # centroid y
    centroidy = wghtTotaly/data.sum()
    centroidydisplay = centroidy +ycoord-radius

    return centroidxdisplay, centroidydisplay


centroidxstar, centroidystar = calculateCentroid(starx, stary, radiusc)
centroidxastrd, centroidyastrd = calculateCentroid(astrdx, astrdy, radiusc)
print centroidxstar, centroidystar
print centroidxastrd, centroidyastrd
deltaxstar = starx - centroidxstar
deltaystar = stary - centroidystar
deltaxastrd = astrdx - centroidxastrd
deltayastrd = astrdy - centroidyastrd
print deltaxstar, deltaystar
print deltaxastrd, deltayastrd

print "star box"
starbox = myimstar[stary-radius-buf:stary+radius+buf+1, starx-radius-buf:starx+radius+buf+1]#extracts a 14 by 14 grid around the given star coordinates
y,x =np.ogrid[-radius-buf:radius+buf+1,-radius-buf:radius+buf+1]
mask = (x+deltaxstar)**2 +(y+deltaystar)**2 >= radius**2
starbox[mask] = 0
print starbox

Nap = np.count_nonzero(starbox)
print Nap

print "total star (star+sky) " # gets the total of all the pixels in the star box
totalstar = np.sum(starbox)
print totalstar








print "blank box"
blankbox = myimblank[stary-radius3-buf:stary+radius3+buf+1, starx-radius3-buf:starx+radius3+buf+1] #extracts a 14 by 14 grid around the given coordinates
y2,x2 =np.ogrid[-radius3-buf:radius3+buf+1,-radius3-buf:radius3+buf+1]
mask2 = (x2+deltaxstar)**2 +(y2+deltaystar)**2 <= (radius2)**2
blankbox[mask2] = 0
mask3 = (x2+deltaxstar)**2 +(y2+deltaystar)**2 >= (radius3)**2
blankbox[mask3]=0
print blankbox

Nbp = np.count_nonzero(blankbox)
print Nbp

print "total blank " # gets the total of all the pixels in the blank box
totalblank = np.sum(blankbox)
print totalblank

print "average blank (avgSky)" # gets the average pixel value in the blank box
averageblank = totalblank/(Nbp)
print averageblank

print "signal star" # calculates the signal for star
signal = totalstar - (averageblank*(Nap))
print signal

print "constant" # gets the value of constant
const = magstar + 2.5*math.log10(signal)
print const







print "asteroid box"
astrdbox = myimastrd[astrdy-radius-buf:astrdy+radius+buf+1, astrdx-radius-buf:astrdx+radius+buf+1] #extracts a 14 by 14 grid around the given asteroid coordinates
y3,x3 =np.ogrid[-radius-buf:radius+buf+1,-radius-buf:radius+buf+1]
mask = (x3+deltaxastrd)**2 +(y3+deltayastrd)**2 >= radius**2
astrdbox[mask] = 0
print astrdbox

print "total astrd (star+sky) " # gets the total of all the pixels in the asteroid box
totalastrd = np.sum(astrdbox)
print totalastrd


print "signal asteroid" # calculates the signal of asteroid
signalastrd = totalastrd - (averageblank*(Nap))
print signalastrd

print "magnitude asteroid" # calculates magnitude of the asteroid
magastrd = -2.5*math.log10(signalastrd) + const
print magastrd





print "instrumental magnitude of sky"
instmag = -2.5*math.log10(averageblank) + const
print instmag

print "sky brightness"
skybright = instmag/1.266**2
print skybright, "mag/arcsec^2"

