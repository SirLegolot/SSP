from __future__ import division
import numpy as np
import matplotlib.pyplot as plt
import math
from astropy.io import fits


#Input image
#im = input ("what is your image? ")
inputimage = fits.getdata("run 3a.fit") # sample image
##myim = fits.getdata("Average1.fit") # jason1984kb

im = inputimage - np.median(inputimage)
#displays image
##plt.imshow(im, vmin=im.mean(), vmax=2*im.mean())
##plt.gray()
##plt.show()

#input coordinates & radius to comput asteroid
##xcoord = input("What is your x coordinate: ")
##ycoord = input("What is your y coordinate: ")
##radius = input("What is your radius to compute centroid: ")
##xcoord = 351
##ycoord = 154
radius = 3


def calculateCentroid(xcoord, ycoord, radius):
    #Creates nine point array centered on the coordinates
    data = im[ycoord-radius:ycoord+radius+1, xcoord-radius:xcoord+radius+1]
##    print data
    

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
##centroidxdisplay, centroidydisplay, errorx, errory = calculateCentroid(xcoord, ycoord, radius)
##print "Centroid (x,y): ", centroidxdisplay, ",", centroidydisplay
##print "Uncertainty(std. dev. of the mean) in x,y: ", errorx, ",", errory

print calculateCentroid(685, 57, 3)
print calculateCentroid(711, 217, 3)
print calculateCentroid(870, 332, 3)
print calculateCentroid(844, 366, 3)
print calculateCentroid(761, 502, 3)
print calculateCentroid(710, 717, 3)
print calculateCentroid(752, 822, 3)
print calculateCentroid(581, 827, 3)
print calculateCentroid(468, 687, 3)
print calculateCentroid(449, 572, 3)
print calculateCentroid(448, 176, 3)
print calculateCentroid(502, 216, 3)
print calculateCentroid(641, 462, 1)

