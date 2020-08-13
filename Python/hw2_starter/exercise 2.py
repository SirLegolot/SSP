from __future__ import division
import numpy as np
import matplotlib.pyplot as plt
import math
from astropy.io import fits

data = np.array([[0,33,21,33,8],
                [0,56,51,53,26],
                [23,120,149,73,18],
                [55,101,116,50,16],
                [11,78,26,2,10]])


# totals list x
totalListx = sum(data)

# weighted total list x
wghtTotalListx = np.copy(totalListx)
for index in range(len(totalListx)):
    wghtTotalListx[index] = wghtTotalListx[index]*index

# weighted total x
wghtTotalx = sum(wghtTotalListx)

# centroid x
centroidx = wghtTotalx/sum(totalListx)
print "centroid x: ", centroidx

# error list x
errorListx = np.copy(totalListx)
for index in range(len(totalListx)):
    errorListx[index] = errorListx[index]*(index-centroidx)**2

# errorx
errorx = math.sqrt(sum(errorListx)/(sum(totalListx)*(sum(totalListx)-1)))
print "error x: ", errorx

# totals list y
totalListy = []
for row in data:
    totalListy.append(sum(row))

# weighted total list y
wghtTotalListy = np.copy(totalListy)
for index in range(len(totalListy)):
    wghtTotalListy[index] = wghtTotalListy[index]*index

# weighted total y
wghtTotaly = sum(wghtTotalListy)

# centroid y
centroidy = wghtTotaly/sum(totalListy)
print "centroid y: ", centroidy

# error list y
errorListy = np.copy(totalListy)
for index in range(len(totalListy)):
    errorListy[index] = errorListy[index]*(index-centroidy)**2

# errorx
errory = math.sqrt(sum(errorListy)/(sum(totalListy)*(sum(totalListy)-1)))
print "error y: ", errory



