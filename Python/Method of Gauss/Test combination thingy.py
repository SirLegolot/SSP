from __future__ import division
import numpy as np
from math import *
import ast

####################################### Initial Setup ###########################################
k=0.01720209895
mu = 1
epsilon = radians(23.4347)
cspeed = 173.145
##NewestDataFile7-19
inputFile = np.genfromtxt("Input3.txt", dtype=None)
##if len(inputFile)>3:
##    print "What combination of observations do you want to compute with?"
##    choices = []
##    for i in range(2, len(inputFile)):
##        print i-1, "-", "[ 1", i , len(inputFile),"]"
##    choice = input("Enter number: ")
##    newinputFile = [inputFile[0], inputFile[choice], inputFile[-1]]
##elif len(inputFile) ==3:
##    newinputFile = inputFile
##print "file\n", newinputFile    


years= []
months = []
days = []
hms = []
RAs = []
Decs = []
Rvectors = []

for item in inputFile:
    item = list(item)
    years.append(item[0])
    months.append(item[1])
    days.append(item[2])
    hms.append([int(number) for number in item[3].split(':')])
    RAs.append(item[4:7])
    Decs.append(item[7:10])
    Rvectors.append(ast.literal_eval(item[10]))

hr = dict()
Min = dict()
sec = dict()
mon = dict()
day = dict()
yr = dict()
RAhr = dict()
RAMin = dict()
RAsec = dict()
deg = dict()
amin = dict()
asec = dict()
for i in range(len(hms)):
    hr[i+1] = hms[i][0]
    Min[i+1] = hms[i][1]
    sec[i+1] = hms[i][2]
for i in range(len(months)):
    mon[i+1] = months[i]
for i in range(len(days)):
    day[i+1] = days[i]
for i in range(len(years)):
    yr[i+1] = years[i]
for i in range(len(RAs)):
    RAhr[i+1] = RAs[i][0]
    RAMin[i+1] = RAs[i][1]
    RAsec[i+1] = RAs[i][2]
for i in range(len(Decs)):
    deg[i+1] = Decs[i][0]
    amin[i+1] = Decs[i][1]
    asec[i+1] = Decs[i][2]


#_______________________Equatorial to Ecliptic______________________________________
def EqtoEc(eqvector):
    matrix = np.array([[1,0,0],
                      [0, cos(epsilon), sin(epsilon)],
                      [0, -sin(epsilon), cos(epsilon)]])
    ecvector = np.dot(matrix, eqvector)
    return ecvector

#_______________________Ecliptic to Equatorial______________________________________
def ECtoEq(ecvector):
    matrix = np.array([[1,0,0],
                      [0, cos(epsilon), -sin(epsilon)],
                      [0, sin(epsilon), cos(epsilon)]])
    eqvector = np.dot(matrix, ecvector)
    return eqvector
Rvector = dict()
for i in range(len(Rvectors)):
    Rvector[i+1] = ECtoEq(np.array(Rvectors[i]))
    Rvector[i+1] = ECtoEq(np.array(Rvectors[i]))
    Rvector[i+1] = ECtoEq(np.array(Rvectors[i]))
print Rvector[1]
