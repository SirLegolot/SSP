from __future__ import division
import numpy as np
from math import *

# 3:47:6.67 6/22/17 UTC --> 21:47:6.67 6/21/17 MDT
# 5:35:34.32 6/26/17 UTC --> 23:35:34.32 6/25/17 MDT
# 5:21:8.219 7/1/17 UTC--> 23:21:8.219 6/30/17 MDT

####################################### Initial Setup ###########################################
k=0.01720209895
mu = 0.01720209895**2
epsilon = radians(23.4347)
cspeed = 173.145

inputFile = np.genfromtxt("OD Data File.txt", dtype=None)

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
    Rvectors.append(item[10:14])
    
hr1 = hms[0][0]
hr2 = hms[1][0]
hr3 = hms[2][0]

Min1 = hms[0][1]
Min2 = hms[1][1]
Min3 = hms[2][1]

sec1 = hms[0][2]
sec2 = hms[1][2]
sec3 = hms[2][2]

mon1 = months[0]
mon2 = months[1]
mon3 = months[2]

day1 = days[0]
day2 = days[1]
day3 = days[2]

yr1 = years[0]
yr2 = years[1]
yr3 = years[2]

RAhr1 = RAs[0][0]
RAhr2 = RAs[1][0]
RAhr3 = RAs[2][0]

RAMin1 = RAs[0][1]
RAMin2 = RAs[1][1]
RAMin3 = RAs[2][1]

RAsec1 = RAs[0][2]
RAsec2 = RAs[1][2]
RAsec3 = RAs[2][2]

deg1 = Decs[0][0]
deg2 = Decs[1][0]
deg3 = Decs[2][0]

amin1 = Decs[0][1]
amin2 = Decs[1][1]
amin3 = Decs[2][1]

asec1 = Decs[0][2]
asec2 = Decs[1][2]
asec3 = Decs[2][2]

Rvector1 = Rvectors[0]
Rvector2 = Rvectors[1]
Rvector3 = Rvectors[2]

#################################### Define general Functions ###################################

#_______________________Find Quadrant____________________________________________________
def findQuadrant(sine, cosine):
    if cosine > 0 and sine > 0: #1
        return asin(sine)

    if cosine < 0 and sine > 0: #2
        return acos(cosine)

    if cosine < 0 and sine < 0: #3
        return pi - asin(sine)

    if cosine > 0 and sine < 0: #4
        return 2*pi + asin(sine)

#_______________________Baby OD____________________________________________________
def babyOD(vecr, vecrdot):
    magvecr = sqrt(vecr[0]**2+vecr[1]**2+vecr[2]**2)

    a = ((2/magvecr)-np.dot(vecrdot, vecrdot))**-1
##    print "a: ",a

    vech = np.cross(vecr, vecrdot)
    magvech = sqrt(vech[0]**2+vech[1]**2+vech[2]**2)
    e = sqrt(1-magvech**2/a)
##    print "e: ", e

    i = acos(vech[2]/magvech)
##    print "i: ", degrees(i)

    sinO = vech[0]/(magvech*sin(i))
    cosO = -vech[1]/(magvech*sin(i))
    O = findQuadrant(sinO, cosO)
##    print "O: ", degrees(O)

    sinU = vecr[2]/(magvecr*sin(i))
    cosU = (vecr[0]*cos(O)+vecr[1]*sin(O))/magvecr
    U = findQuadrant(sinU, cosU)

    sinV = (a*(1-e**2)*np.dot(vecr, vecrdot))/(magvech*magvecr*e)
    cosV = ((a*(1-e**2)/magvecr)-1)/e
    V = findQuadrant(sinV, cosV)
##    print "V", V
    print "cos(E)", 1-magvecr/a/e
    w = U-V
    if w<0:
        w= w+2*pi
##    print magvecr
    
    E = acos((1-magvecr/a)*(1/e))
    if V<0:
        E = E- 2*pi
    elif V>2*pi:
        E = E- 2*pi
    elif 0<V<pi:
        E=E
    elif pi< V < 2*pi:
        E = 2*pi-E
        
    M = E - e*sin(E)
    

    return a, e, i, O, w, M, E

#_______________________Equatorial to Ecliptic______________________________________
def EqtoEc(eqvector):
    matrix = np.array([[1,0,0],
                      [0, cos(epsilon), sin(epsilon)],
                      [0, -sin(epsilon), cos(epsilon)]])
    ecvector = np.dot(matrix, eqvector)
    return ecvector

#_______________________Find E from M function______________________________________
def solvekep(Mtrue):
    E = Mtrue 
    Mguess = E - e*sin(E)
    while abs(Mguess - Mtrue) > 10**-4: 
        Mguess = E - e*sin(E)
        E = E - (Mguess - Mtrue) / (1 - e*cos(E))
    return E

#_______________________Convert Civil Date to Julian Date___________________________

def CDtoJD(hr, Min, sec, yr, mon, day):
    J0 = 367*yr - int((7*(yr+int((mon+9)/12)))/4) + int(275*mon/9) + day + 1721013.5
    JD = J0 + (hr+Min/60+sec/3600)/24
    return JD

#_______________________Convert RA to radians_______________________________________
def radRA(hr, Min, sec):
    RArad = radians((hr+Min/60+sec/3600)*15)
    return RArad

#_______________________Convert Dec to radians______________________________________
def radDec(deg, amin, asec):
    Decrad = radians(deg + amin/60 +asec/3600)
    if deg<0:
        Decrad = radians(deg - amin/60 -asec/3600)
    return Decrad

################################### First Iteration #############################################

#_______________________Initial Time Stuff__________________________________________
t1 = CDtoJD(hr1, Min1, sec1, yr1, mon1, day1)
t2 = CDtoJD(hr2, Min2, sec2, yr2, mon2, day2)
t3 = CDtoJD(hr3, Min3, sec3, yr3, mon3, day3)
tau3 = k*(t3-t2)
tau1 = k*(t1-t2)
tau = tau3 - tau1
##print t1, t2, t3
##print tau3, tau1, tau

#_______________________Convert RA to radians_______________________________________
RA1 = radRA(RAhr1, RAMin1, RAsec1)
RA2 = radRA(RAhr2, RAMin2, RAsec2)
RA3 = radRA(RAhr3, RAMin3, RAsec3)
##print "RAs:", RA1, RA2 ,RA3

#_______________________Convert Dec to radians______________________________________
Dec1 = radDec(deg1, amin1, asec1)
Dec2 = radDec(deg2, amin2, asec2)
Dec3 = radDec(deg3, amin3, asec3)
##print "Dec:", Dec1, Dec2, Dec3

#_______________________Convert RA and Dec to Unit Rho______________________________
def ADtoUnitRho(RA, Dec):
    unitRho = np.array((cos(RA)*cos(Dec), sin(RA)*cos(Dec), sin(Dec)))
    return unitRho
unitRho1 = ADtoUnitRho(RA1, Dec1)
unitRho2 = ADtoUnitRho(RA2, Dec2)
unitRho3 = ADtoUnitRho(RA3, Dec3)
##unitRho1 = EqtoEc(unitRho1)
##unitRho2 = EqtoEc(unitRho2)
##unitRho3 = EqtoEc(unitRho3)
##print unitRho1, unitRho2, unitRho3

#_______________________D stuff_____________________________________________________
D11 = np.dot(np.cross(Rvector1, unitRho2), unitRho3)
##print "D11", D11
D12 = np.dot(np.cross(Rvector2, unitRho2), unitRho3)
##print "D12", D12
D13 = np.dot(np.cross(Rvector3, unitRho2), unitRho3)
##print "D13", D13
D21 = np.dot(np.cross(unitRho1, Rvector1), unitRho3)
##print "D21", D21
D22 = np.dot(np.cross(unitRho1, Rvector2), unitRho3)
##print "D22", D22
D23 = np.dot(np.cross(unitRho1, Rvector3), unitRho3)
##print "D23", D23
D31 = np.dot(unitRho1, np.cross(unitRho2,Rvector1))
##print "D31", D31
D32 = np.dot(unitRho1, np.cross(unitRho2,Rvector2))
##print "D32", D32
D33 = np.dot(unitRho1, np.cross(unitRho2,Rvector3))
##print "D33", D33
D0 = np.dot(unitRho1, np.cross(unitRho2, unitRho3))
##print "D0",D0

#_______________________Intermediate Steps to get R2________________________________
A1 = tau3/tau
A3 = -tau1/tau
A = -(A1*D21 - D22 + A3*D23)/D0
B1 = (A1/6)*(tau**2-tau3**2)
B3 = (A3/6)*(tau**2-tau1**2)
B = -(B1*D21 + B3*D23)/D0
E = -2*(np.dot(unitRho2, Rvector2))
F = (sqrt(Rvector2[0]**2 + Rvector2[1]**2 + Rvector2[2]**2))**2
##print A, B, E, F
##print A1, B1, A3, B3
a = -(A**2 + A*E + F)##print tau3, tau1, tau

b = -mu*(2*A*B + B*E)
c = -mu**2 * B**2
##print a, b, c

#_______________________Solve equation to find r2__________________________________
roots = np.roots([1, 0, a, 0, 0, b, 0, 0, c])
##print roots
realRootsbool = np.isreal(roots)
##print realRootsbool
realRoots = []
for i in range(len(realRootsbool)):
    if realRootsbool[i] == True:
        realRoots.append(roots[i])
realRoots = np.real(realRoots) # there could be multiple roots
##print realRoots##print tau3, tau1, tau

possibleRoots = []
for i in range(len(realRoots)):
    if realRoots[i] > 0.1 and realRoots[i]<20:
        possibleRoots.append(realRoots[i])
r2 = possibleRoots[0] # must select a root
##print "possibleRoots:", possibleRoots
##print "r2:",r2

#_______________________Initial fi and gi___________________________________________
f1 = 1-((mu*tau1)/(2*r2**3))
g1 = tau1-((mu*tau1**3)/(6*r2**3))
f3 = 1-((mu*tau3)/(2*r2**3))
g3 = tau3-((mu*tau3**3)/(6*r2**3))
print f1, g1, f3, g3

#_______________________Get rho mag values__________________________________________
def getRhoMags():
    C1 = g3/(f1*g3-g1*f3)
    C2=-1
    C3 = -g1/(f1*g3-g1*f3)
##    print C1, C3
##    print D11, D21, D31, D12, D22, D23, D13, D23, D33, D0
    rho1mag = (C1*D11+C2*D12+C3*D13)/(C1*D0)
    rho2mag = (C1*D21+C2*D22+C3*D23)/(C2*D0)
    rho3mag = (C1*D31+C2*D32+C3*D33)/(C3*D0)
    print rho1mag, rho2mag, rho3mag
    return rho1mag, rho2mag, rho3mag
rho1mag, rho2mag, rho3mag = getRhoMags()
##print rho1mag, rho2mag, rho3mag

#_______________________Get r vectors and r dot vectors_____________________________
def getrvectors():
    rvector1 = (rho1mag*unitRho1)-Rvector1
    rvector2 = (rho2mag*unitRho2)-Rvector2
    rvector3 = (rho3mag*unitRho3)-Rvector3
##    print rvector2
    d1 = -f3/(f1*g3-f3*g1)
    d3 = f1/(f1*g3 - f3*g1)
    rDotVector2 = (d1*rvector1)+(d3*rvector3)
    return rvector1, rvector2, rvector3, rDotVector2
rvector1, rvector2, rvector3, rDotVector2 = getrvectors()
##print "rvector2", rvector2
##print "rDotVector2:", rDotVector2

#_______________________Find elements for r2v and r2dotv____________________________
a2, e2, i2, O2, w2, M2, E2 = babyOD(rvector2, rDotVector2)
##print "a2:", a2
##print "e2:", e2
##print "i2:", i2
##print "O2:", O2
##print "w2:", w2
##print "M2:", M2
##print "E2:", E2

################################### Second Iteration ##########################################
###_______________________Light Corrections_______________________________
##t1 = t1-(rho1mag/cspeed)
##t2 = t2-(rho2mag/cspeed)
##t3 = t3-(rho3mag/cspeed)
##tau3 = k*(t3-t2)
##tau1 = k*(t1-t2)
##tau = tau3 - tau1
####    print tau3, tau1, tau
####    print t1, t2, t3
##
###_______________________Calculate E1, E3 using newton's method__________
##n = a2**(-3/2)
##M1 = M2 + n*(t1-t2)
##M3 = M2 + n*(t3-t2)
####    print M1, M3
##E1 = solvekep(M1)
##E3 = solvekep(M3)
####    print E1, E3
##
###_______________________Calculate new fi and gi_________________________
####f1 = 1-a2/r2*(1-cos(E1-E2))
####f3 = 1-a2/r2*(1-cos(E3-E2))
####g1 = tau1+(1/n)*(sin(E1-E2)-(E1-E2))
####g3 = tau3+(1/n)*(sin(E3-E2)-(E3-E2))
##f1 = 1-((mu*tau1)/(2*r2**3))+mu*(np.dot(rvector2, rDotVector2))*tau1**3/(2*r2**5)
##g1 = tau1-((mu*tau1**3)/(6*r2**3))
##f3 = 1-((mu*tau3)/(2*r2**3))+mu*(np.dot(rvector2, rDotVector2))*tau3**3/(2*r2**5)
##g3 = tau3-((mu*tau3**3)/(6*r2**3))
####    print f1, g1, f3, g3
##(R
###_______________________Get new rho values______________________________
##rho1mag, rho2mag, rho3mag = getRhoMags()
####    print rho1mag, rho2mag, rho3mag
##
##
###_______________________Get new r vectors and r dot vectors_____________
##rvector1, rvector2, rvector3, rDotVector2 = getrvectors()
##r2 = np.linalg.norm(rvector2)
##
###_______________________Baby OD_________________________________________
##a2, e2, i2, O2, w2, M2, E2 = babyOD(rvector2, rDotVector2)
################################## Loop Iteration ###############################################
counter = 0
while True:
    #_______________________Light Corrections_______________________________
    t1 = t1-(rho1mag/cspeed)
    t2 = t2-(rho2mag/cspeed)
    t3 = t3-(rho3mag/cspeed)
    tau3 = k*(t3-t2)
    tau1 = k*(t1-t2)
    tau = tau3 - tau1
##    print tau3, tau1, tau
##    print t1, t2, t3

    #_______________________Calculate E1, E3 using newton's method__________
    n = a2**(-3/2)
    M1 = M2 + n*(t1-t2)
    M3 = M2 + n*(t3-t2)
##    print M1, M3
    E1 = solvekep(M1)
    E3 = solvekep(M3)
##    print E1, E3

    #_______________________Calculate new fi and gi_________________________
    f1 = 1-a2/r2*(1-cos(E1-E2))
    f3 = 1-a2/r2*(1-cos(E3-E2))
    g1 = tau1+(1/n)*(sin(E1-E2)-(E1-E2))
    g3 = tau3+(1/n)*(sin(E3-E2)-(E3-E2))
####    f1 = 1-((mu*tau1)/(2*r2**3))+mu*(np.dot(rvector2, rDotVector2))*tau1**3/(2*r2**5)
####    g1 = tau1-((mu*tau1**3)/(6*r2**3))
####    f3 = 1-((mu*tau3)/(2*r2**3))+mu*(np.dot(rvector2, rDotVector2))*tau3**3/(2*r2**5)
####    g3 = tau3-((mu*tau3**3)/(6*r2**3))
##    u = mu/r2**3
##    z = np.dot(rvector2, rDotVector2)/r2**2
##    q = np.dot(rDotVector2, rDotVector2)/r2**2 - u
##    f1 = 1-((mu*tau1)/(2*r2**3))+mu*(np.dot(rvector2, rDotVector2))*tau1**3/(2*r2**5) +(3*u*q - 15*u*z**2 + u**2)*(tau1**4)/24
##    g1 = tau1-((mu*tau1**3)/(6*r2**3)) + u*z*(tau1**4)/4
##    f3 = 1-((mu*tau3)/(2*r2**3))+mu*(np.dot(rvector2, rDotVector2))*tau3**3/(2*r2**5) +(3*u*q - 15*u*z**2 + u**2)*(tau3**4)/24
##    g3 = tau3-((mu*tau3**3)/(6*r2**3)) + u*z*(tau3**4)/4
    print f1, g1, f3, g3

    #_______________________Get new rho values______________________________
    rho1mag, rho2mag, rho3mag = getRhoMags()
##    print rho1mag, rho2mag, rho3mag
 

    #_______________________Get new r vectors and r dot vectors_____________
    rvector1, rvector2, rvector3, rDotVector2 = getrvectors()
    r2old = r2
    r2 = np.linalg.norm(rvector2)

    #_______________________Baby OD_________________________________________
    a2, e2, i2, O2, w2, M2, E2 = babyOD(rvector2, rDotVector2)
    counter +=1

    #_______________________Condition of breaking loop______________________
    if abs(r2-r2old)<10**-12:
        break


    
print "new r2", r2
print "new rvector2", rvector2
print "new rDotVector2:", rDotVector2




