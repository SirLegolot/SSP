from __future__ import division
import numpy as np
import matplotlib.pyplot as plt
import math
from astropy.io import fits


#Input image
#im = input ("what is your image? ")
im = fits.getdata("sampleimage.fits") # sample image
myim = fits.getdata("Average1.fit") # jason1984kb

#displays image
plt.imshow(im, vmin=im.mean(), vmax=2*im.mean())
plt.gray()
plt.show()

#input coordinates & radius to comput asteroid
##xcoord = input("What is your x coordinate: ")
##ycoord = input("What is your y coordinate: ")
##radius = input("What is your radius to compute centroid: ")
xcoord = 351
ycoord = 154
radius = 1


def calculateCentroid(xcoord, ycoord, radius):
    #Creates nine point array centered on the coordinates
    data = im[ycoord-radius:ycoord+radius+1, xcoord-radius:xcoord+radius+1]
    

    # totals list x
    totalListx = 0
    for i in data[:,:]:
        totalListx = totalListx + i
    print totalListx

    # weighted total x
    wghtTotalx = 0
    for index in range(len(totalListx)):
        wghtTotalx = totalListx[index]*index + wghtTotalx

    # centroid x
    centroidx = wghtTotalx/data.sum()
    centroidxdisplay = centroidx +xcoord-1

    # error list x
    errorsumx = 0
    for index in range(len(totalListx)):
        errorsumx = totalListx[index]*(index-centroidx)**2 + errorsumx

    # errorx
    errorx = math.sqrt(errorsumx/(data.sum()*(data.sum()-1)))

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
    centroidydisplay = centroidy +ycoord-1

    # error list y
    errorsumy = 0
    for index in range(len(totalListy)):
        errorsumy = totalListy[index]*(index-centroidy)**2 + errorsumy

    # errorx
    errory = math.sqrt(errorsumy/(data.sum()*(data.sum()-1)))
    return centroidxdisplay, centroidydisplay, errorx, errory

centroidxdisplay, centroidydisplay, errorx, errory = calculateCentroid(xcoord, ycoord, radius)
print "Centroid (x,y): ", centroidxdisplay, ",", centroidydisplay
print "Uncertainty(std. dev. of the mean) in x,y: ", errorx, ",", errory

