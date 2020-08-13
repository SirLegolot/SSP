from __future__ import division
from math import *
from astropy.io import fits
import numpy as np

dark1 = fits.getdata('cloudyDark00000001.fit')[20:-20,20:-20]
dark2 = fits.getdata('cloudyDark00000003.fit')[20:-20,20:-20]
dark3 = fits.getdata('cloudyDark00000005.fit')[20:-20,20:-20]
dark4 = fits.getdata('cloudyDark00000007.fit')[20:-20,20:-20]
dark5 = fits.getdata('cloudyDark00000009.fit')[20:-20,20:-20]
dark6 = fits.getdata('cloudyDark00000011.fit')[20:-20,20:-20]
dark7 = fits.getdata('cloudyDark00000013.fit')[20:-20,20:-20]

bias1 = fits.getdata('cloudyBias00000002.fit')[20:-20,20:-20]
bias2 = fits.getdata('cloudyBias00000004.fit')[20:-20,20:-20]
bias3 = fits.getdata('cloudyBias00000006.fit')[20:-20,20:-20]
bias4 = fits.getdata('cloudyBias00000008.fit')[20:-20,20:-20]
bias5 = fits.getdata('cloudyBias00000010.fit')[20:-20,20:-20]
bias6 = fits.getdata('cloudyBias00000012.fit')[20:-20,20:-20]
bias7 = fits.getdata('cloudyBias00000014.fit')[20:-20,20:-20]

#darks - bias
img1 = dark1 - bias1
img2 = dark2 - bias2
img3 = dark3 - bias3
img4 = dark4 - bias4
img5 = dark5 - bias5
img6 = dark6 - bias6
img7 = dark7 - bias7

#signal for corrected dark frames
print "signal for corrected dark frames as temp changes"
##signal1 = np.average(img1)
##signal2 = np.average(img2)
##signal3 = np.average(img3)
##signal4 = np.average(img4)
##signal5 = np.average(img5)
##signal6 = np.average(img6)
##signal7 = np.average(img7)
##print signal1
##print signal2
##print signal3
##print signal4
##print signal5
##print signal6
##print signal7



#_________________________________Bias stuff_____________________________
#average bias
avgbias1 = np.average(bias1)
avgbias2 = np.average(bias2)
avgbias3 = np.average(bias3)
avgbias4 = np.average(bias4)
avgbias5 = np.average(bias5)
avgbias6 = np.average(bias6)
avgbias7 = np.average(bias7)

##print "\nsignal for bias frames over time as temp changes"
##print avgbias1
##print avgbias2
##print avgbias3
##print avgbias4
##print avgbias5
##print avgbias6
##print avgbias7

#standard deviation of the pixel counts for bias
sigmab1 =  np.std(bias1)
sigmab2 =  np.std(bias2)
sigmab3 =  np.std(bias3)
sigmab4 =  np.std(bias4)
sigmab5 =  np.std(bias5)
sigmab6 =  np.std(bias6)
sigmab7 =  np.std(bias7)




#dark count rate
print "\ndark count rate as temp changes"
ratedark1 = signal1/600
ratedark2 = signal2/600
ratedark3 = signal3/600
ratedark4 = signal4/600
ratedark5 = signal5/600
ratedark6 = signal6/600
ratedark7 = signal7/600

print ratedark1
print ratedark2
print ratedark3
print ratedark4
print ratedark5
print ratedark6
print ratedark7




###avg dark
##avgdark1 = np.average(dark1)
##avgdark2 = np.average(dark2)
##avgdark3 = np.average(dark3)
##avgdark4 = np.average(dark4)
##avgdark5 = np.average(dark5)
##avgdark6 = np.average(dark6)
##avgdark7 = np.average(dark7)
##
##print "\navg values for dark frames"
##print avgdark1
##print avgdark2
##print avgdark3
##print avgdark4
##print avgdark5
##print avgdark6
##print avgdark7


###Standard Error propagation for bias subtracted frames
##sigma12 = np.average(abs((img1 - img2))/sqrt(2))
##sigma23 = np.average(abs((img2 - img3))/sqrt(2))
##sigma34 = np.average(abs((img3 - img4))/sqrt(2))
##sigma45 = np.average(abs((img4 - img5))/sqrt(2))
##sigma56 = np.average(abs((img5 - img6))/sqrt(2))
##sigma67 = np.average(abs((img6 - img7))/sqrt(2))
##print "\nnoise levels for bias subtracted frames:"
##print sigma12
##print sigma23
##print sigma34
##print sigma45
##print sigma56
##print sigma67


###Standard Error propagation for bias
##b_sigma12 = np.average(abs((bias1 - bias2))/sqrt(2))
##b_sigma23 = np.average(abs((bias2 - bias3))/sqrt(2))
##b_sigma34 = np.average(abs((bias3 - bias4))/sqrt(2))
##b_sigma45 = np.average(abs((bias4 - bias5))/sqrt(2))
##b_sigma56 = np.average(abs((bias5 - bias6))/sqrt(2))
##b_sigma67 = np.average(abs((bias6 - bias7))/sqrt(2))

##print "\nnoise levels for bias frames:"
##print b_sigma12
##print b_sigma23
##print b_sigma34
##print b_sigma45
##print b_sigma56
##print b_sigma67

##d_sigma12 = np.average(abs((dark1 - dark2))/sqrt(2))
##d_sigma23 = np.average(abs((dark2 - dark3))/sqrt(2))
##d_sigma34 = np.average(abs((dark3 - dark4))/sqrt(2))
##d_sigma45 = np.average(abs((dark4 - dark5))/sqrt(2))
##d_sigma56 = np.average(abs((dark5 - dark6))/sqrt(2))
##d_sigma67 = np.average(abs((dark6 - dark7))/sqrt(2))
##
##print "\nnoise levels for dark frames:"
##print d_sigma12
##print d_sigma23
##print d_sigma34
##print d_sigma45
##print d_sigma56
##print d_sigma67



