from __future__ import division
from math import *
from astropy.io import fits
import numpy as np


bias1 = fits.getdata('cloudyBias00000002.fit')[20:-20,20:-20]
bias2 = fits.getdata('cloudyBias00000004.fit')[20:-20,20:-20]
bias3 = fits.getdata('cloudyBias00000006.fit')[20:-20,20:-20]
bias4 = fits.getdata('cloudyBias00000008.fit')[20:-20,20:-20]
bias5 = fits.getdata('cloudyBias00000010.fit')[20:-20,20:-20]
bias6 = fits.getdata('cloudyBias00000012.fit')[20:-20,20:-20]
bias7 = fits.getdata('cloudyBias00000014.fit')[20:-20,20:-20]


#average bias
avgbias1 = np.average(bias1)
avgbias2 = np.average(bias2)
avgbias3 = np.average(bias3)
avgbias4 = np.average(bias4)
avgbias5 = np.average(bias5)
avgbias6 = np.average(bias6)
avgbias7 = np.average(bias7)

#/sqrt(2*984**2-2)
#standard deviation of the pixel counts for bias
##sigmab1 =  np.std(bias1)
##sigmab2 =  np.std(bias2)
##sigmab3 =  np.std(bias3)
##sigmab4 =  np.std(bias4)
##sigmab5 =  np.std(bias5)
##sigmab6 =  np.std(bias6)
##sigmab7 =  np.std(bias7)

signals = [avgbias1, avgbias2, avgbias3, avgbias4, avgbias5, avgbias6, avgbias7]
print signals

sigmab1 =  np.std(signals)
sigmab2 =  np.std(signals)
sigmab3 =  np.std(signals)
sigmab4 =  np.std(signals)
sigmab5 =  np.std(signals)
sigmab6 =  np.std(signals)
sigmab7 =  np.std(signals)
##print sigmab1

temps = [292.5182-273, 288.4548-273, 285.0937-273, 280.9744-273, 276.5926-273, 272.7682-273, 269.0117-273]
sigmas = [sigmab1, sigmab2, sigmab3, sigmab4, sigmab5, sigmab6, sigmab7]

#temp = x, signal = y
oneOverSigmaB2 = 0
for i in range(len(sigmas)):
    oneOverSigmaB2 = oneOverSigmaB2 + (1/sigmas[i])**2
x2OverSigmaB2 = 0
for i in range(len(sigmas)):
    x2OverSigmaB2 = x2OverSigmaB2 + (temps[i]/sigmas[i])**2
xOverSigmaB2 = 0
for i in range(len(sigmas)):
    xOverSigmaB2 = xOverSigmaB2 + temps[i]/(sigmas[i])**2
yOverSigmaB2 = 0
for i in range(len(sigmas)):
    yOverSigmaB2 = yOverSigmaB2 + signals[i]/(sigmas[i])**2
xyOverSigmaB2 = 0
for i in range(len(sigmas)):
    xyOverSigmaB2 = xyOverSigmaB2 + temps[i]*signals[i]/(sigmas[i])**2

Delta = oneOverSigmaB2*x2OverSigmaB2-xOverSigmaB2**2
a = (1/Delta)*(xOverSigmaB2*yOverSigmaB2-xOverSigmaB2*xyOverSigmaB2)
b = (1/Delta)*(oneOverSigmaB2*xyOverSigmaB2-xOverSigmaB2*yOverSigmaB2)
print Delta, a, b

    
    
