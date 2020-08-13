from __future__ import division
import numpy as np
from math import *
import ast
from Ephemerisgeneration import ephem

####################################### Initial Setup ###########################################
k=0.01720209895
mu = 1
epsilon = radians(23.4347)
cspeed = 173.145
##NewestDataFile7-19 AlexTest
inputFile = np.genfromtxt("PradeepInput.txt", dtype=None)
##print inputFile

# This creates a new list with the data from the three observations that are selected by the user
if len(inputFile)>3:
    print "What combination of observations do you want to compute with?"
    for i in range(2, len(inputFile)):
        print "("+str(i-1)+")", "[ 1", i , len(inputFile),"]"
    choice = input("Enter number: ")
    newinputFile = [inputFile[0], inputFile[choice], inputFile[-1]]
elif len(inputFile) ==3:
    newinputFile = inputFile
##print newinputFile
# 1= closed form functions, 2 = 3rd order taylor series, 3 = 4th order taylor series
fg = 1
    
    
years= []
months = []
days = []
hms = []
RAs = []
Decs = []
Rvectors = []

for item in newinputFile:
    item = list(item)
    years.append(item[0])
    months.append(item[1])
    days.append(item[2])
    hms.append([int(number) for number in item[3].split(':')])
    RAs.append(item[4:7])
    Decs.append(item[7:10])
##    Rvectors.append(item[10:13])
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
##print Rvector[1], Rvector[2], Rvector[3]
#################################### Define general Functions ###################################

#_______________________Find Quadrant_______________________________________________
def findQuadrant(sine, cosine):
    if cosine > 0 and sine > 0: #1
        return asin(sine)

    if cosine < 0 and sine > 0: #2
        return acos(cosine)

    if cosine < 0 and sine < 0: #3
        return pi - asin(sine)

    if cosine > 0 and sine < 0: #4
        return 2*pi + asin(sine)

#_______________________Get a_______________________________________________________
def geta(vecr, vecrdot):
    magvecr = np.linalg.norm(vecr)
##    print magvecr
    a = ((2/magvecr)-np.dot(vecrdot, vecrdot))**-1
##    print "a: ",a
    return a
#_______________________Baby OD____________________________________________________
def babyOD(vecr, vecrdot):
    magvecr = np.linalg.norm(vecr)

    a = ((2/magvecr)-np.dot(vecrdot, vecrdot))**-1
##    print "a: ",a

    vech = np.cross(vecr, vecrdot)
    magvech = np.linalg.norm(vech)
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
##    print "V: ", V
##    print "cos(E)", (1-magvecr/a)*(1/e)
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
    M  = M + (2457956.75-t2)*k*n   # use 2456453.5 for jpl comparison, currently for july 22, 2017 6:00 UT
    M = M % (2*pi)

    

    return a, e, i, O, w, M, E

#_______________________Newton's Method_____________________________________________
def newton(tau_i):
    Enot = n*tau_i
    E =  n*tau_i
    f = E - (1-r2/a2)*sin(E) + np.dot(rvector2,rDotVector2)*(1-cos(E))/(n*a2**2)-Enot
    fprime = 1- (1-r2/a2)*cos(E) + np.dot(rvector2,rDotVector2)*sin(E)/(n*a2**2)
    while abs(f/fprime) > 10**-4: 
        E = E - f/fprime
        f = E - (1-r2/a2)*sin(E) + np.dot(rvector2,rDotVector2)*(1-cos(E))/(n*a2**2)-Enot
        fprime = 1- (1-r2/a2)*cos(E) + np.dot(rvector2,rDotVector2)*sin(E)/(n*a2**2)
##        print E
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
t1orig = CDtoJD(hr[1], Min[1], sec[1], yr[1], mon[1], day[1])
t2orig = CDtoJD(hr[2], Min[2], sec[2], yr[2], mon[2], day[2])
t3orig = CDtoJD(hr[3], Min[3], sec[3], yr[3], mon[3], day[3])
t1 = CDtoJD(hr[1], Min[1], sec[1], yr[1], mon[1], day[1])
t2 = CDtoJD(hr[2], Min[2], sec[2], yr[2], mon[2], day[2])
t3 = CDtoJD(hr[3], Min[3], sec[3], yr[3], mon[3], day[3])
tau3 = k*(t3-t2)
tau1 = k*(t1-t2)
tau = tau3 - tau1
##print t1, t2, t3
##print tau3, tau1, tau

#_______________________Convert RA to radians_______________________________________
RA1 = radRA(RAhr[1], RAMin[1], RAsec[1])
RA2 = radRA(RAhr[2], RAMin[2], RAsec[2])
RA3 = radRA(RAhr[3], RAMin[3], RAsec[3])
##print "RAs:", RA1, RA2 ,RA3

#_______________________Convert Dec to radians______________________________________
Dec1 = radDec(deg[1], amin[1], asec[1])
Dec2 = radDec(deg[2], amin[2], asec[2])
Dec3 = radDec(deg[3], amin[3], asec[3])
##print "Dec:", Dec1, Dec2, Dec3

#_______________________Convert RA and Dec to Rho hats______________________________
def ADtoUnitRho(RA, Dec):
    unitRho = np.array((cos(RA)*cos(Dec), sin(RA)*cos(Dec), sin(Dec)))
    return unitRho
unitRho1 = ADtoUnitRho(RA1, Dec1)
unitRho2 = ADtoUnitRho(RA2, Dec2)
unitRho3 = ADtoUnitRho(RA3, Dec3)
##print unitRho1, unitRho2, unitRho3

#_______________________D stuff_____________________________________________________
D11 = np.dot(np.cross(Rvector[1], unitRho2), unitRho3)
D12 = np.dot(np.cross(Rvector[2], unitRho2), unitRho3)
D13 = np.dot(np.cross(Rvector[3], unitRho2), unitRho3)
D21 = np.dot(np.cross(unitRho1, Rvector[1]), unitRho3)
D22 = np.dot(np.cross(unitRho1, Rvector[2]), unitRho3)
D23 = np.dot(np.cross(unitRho1, Rvector[3]), unitRho3)
D31 = np.dot(unitRho1, np.cross(unitRho2,Rvector[1]))
D32 = np.dot(unitRho1, np.cross(unitRho2,Rvector[2]))
D33 = np.dot(unitRho1, np.cross(unitRho2,Rvector[3]))
D0 = np.dot(unitRho1, np.cross(unitRho2, unitRho3))
##print "D0 ", D0
##print "D11 ", D11
##print "D12 ", D12
##print "D13 ", D13
##print "D21 ", D21
##print "D22 ", D22
##print "D23 ", D23
##print "D31 ", D31
##print "D32 ", D32
##print "D33 ", D33


#_______________________Intermediate Steps to get R2________________________________
A1 = tau3/tau
A3 = -tau1/tau
A = -(A1*D21 - D22 + A3*D23)/D0
B1 = (A1/6)*(tau**2-tau3**2)
B3 = (A3/6)*(tau**2-tau1**2)
B = -(B1*D21 + B3*D23)/D0
E = -2*(np.dot(unitRho2, Rvector[2]))
F = (np.linalg.norm(Rvector[2]))**2
##print A, B, E, F
##print A1, B1, A3, B3
a = -(A**2 + A*E + F)
b = -mu*(2*A*B + B*E)
c = -mu**2 * B**2
##print a, b, c

#_______________________Solve equation to find r2__________________________________
roots = np.roots([1, 0, a, 0, 0, b, 0, 0, c])
##print "roots:", roots

realRootsbool = np.isreal(roots)
##print "realRootsbool: ", realRootsbool
realRoots = []
for i in range(len(realRootsbool)):
    if realRootsbool[i] == True:
        realRoots.append(roots[i])
realRoots = np.real(realRoots) # there could be multiple roots
##print "real Roots: ", realRoots

possibleRoots = []
for i in range(len(realRoots)):
    if realRoots[i] > 0.1 and realRoots[i]<20:
        possibleRoots.append(realRoots[i])
if len(possibleRoots)>1:
    print "\nPossible Roots:"
    for i in range(len(possibleRoots)):
        print "("+str(i+1)+")",possibleRoots[i]
    choose = input("Which root do you want to choose: ")
    r2 = possibleRoots[choose-1]
else:
    r2 = possibleRoots[0]
##print "initial r2:",r2

#_______________________Initial fi and gi___________________________________________
f1 = 1-((mu*tau1**2)/(2*r2**3))
g1 = tau1-((mu*tau1**3)/(6*r2**3))
f3 = 1-((mu*tau3**2)/(2*r2**3))
g3 = tau3-((mu*tau3**3)/(6*r2**3))
##print f1, g1, f3, g3

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
##    print "Rho mags:", rho1mag, rho2mag, rho3mag
    return rho1mag, rho2mag, rho3mag
rho1mag, rho2mag, rho3mag = getRhoMags()
##print rho1mag, rho2mag, rho3mag

#_______________________Get r vectors and r dot vectors_____________________________
def getrvectors():
    rvector1 = (rho1mag*unitRho1)-Rvector[1]
    rvector2 = (rho2mag*unitRho2)-Rvector[2]
    rvector3 = (rho3mag*unitRho3)-Rvector[3]
    r2 = np.linalg.norm(rvector2)
##    print r2
    d1 = -f3/(f1*g3-f3*g1)
    d3 = f1/(f1*g3 - f3*g1)
    rDotVector2 = (d1*rvector1)+(d3*rvector3)
    r2dot = np.linalg.norm(rDotVector2)
##    print r2dot
    return rvector1, rvector2, rvector3, rDotVector2, r2, r2dot
rvector1, rvector2, rvector3, rDotVector2, r2, r2dot = getrvectors()
##print "intital rvector2", rvector2
##print "intital rDotVector2:", rDotVector2

#_______________________Find elements for r2v and r2dotv____________________________
a2= geta(rvector2, rDotVector2)
##print "initial a2:", a2


################################## Loop Iteration ###############################################
converge = 0
counter = 0
while True:
    #_______________________Light Corrections_______________________________
    t1 = t1orig-(rho1mag/cspeed)
    t2 = t2orig-(rho2mag/cspeed)
    t3 = t3orig-(rho3mag/cspeed)
    tau3 = k*(t3-t2)
    tau1 = k*(t1-t2)
    tau = tau3 - tau1
##    print tau3, tau1, tau
##    print t1, t2, t3

    #_______________________Calculate new fs and gs_________________________
    n = a2**(-3/2)
    if fg == 1:
        deltaE1 = newton(tau1)
        deltaE3 = newton(tau3)
        f1 = 1-a2/r2*(1-cos(deltaE1))                                             
        f3 = 1-a2/r2*(1-cos(deltaE3))                                           
        g1 = tau1+(1/n)*(sin(deltaE1)-(deltaE1))
        g3 = tau3+(1/n)*(sin(deltaE3)-(deltaE3))
    elif fg ==2:
        f1 = 1-((mu*tau1**2)/(2*r2**3))+mu*(np.dot(rvector2, rDotVector2))*tau1**3/(2*r2**5)
        g1 = tau1-((mu*tau1**3)/(6*r2**3))
        f3 = 1-((mu*tau3**2)/(2*r2**3))+mu*(np.dot(rvector2, rDotVector2))*tau3**3/(2*r2**5)
        g3 = tau3-((mu*tau3**3)/(6*r2**3))
    elif fg ==3:
        u = mu/r2**3
        z = np.dot(rvector2, rDotVector2)/r2**2
        q = np.dot(rDotVector2, rDotVector2)/r2**2 - u
        f1 = 1-((mu*tau1**2)/(2*r2**3))+mu*(np.dot(rvector2, rDotVector2))*tau1**3/(2*r2**5) +(3*u*q - 15*u*z**2 + u**2)*(tau1**4)/24
        g1 = tau1-((mu*tau1**3)/(6*r2**3)) + u*z*(tau1**4)/4
        f3 = 1-((mu*tau3**2)/(2*r2**3))+mu*(np.dot(rvector2, rDotVector2))*tau3**3/(2*r2**5) +(3*u*q - 15*u*z**2 + u**2)*(tau3**4)/24
        g3 = tau3-((mu*tau3**3)/(6*r2**3)) + u*z*(tau3**4)/4
##    print E1, E3
##    print f1, g1, f3, g3

    #_______________________Get new rho values______________________________
    rho1mag, rho2mag, rho3mag = getRhoMags()
##    print rho1mag, rho2mag, rho3mag

    #_______________________Get new r vectors and r dot vectors_____________
    r2old = r2
    rvector1, rvector2, rvector3, rDotVector2, r2, r2dot = getrvectors()
##    print "r2",r2
    

    #_______________________Get new a_________________________________________
    a2= geta(rvector2, rDotVector2)

    #_______________________Condition of breaking loop______________________
    if abs(r2-r2old)<10**-12:
        converge = 1
        break
    elif counter>1000:
        print "Did not converge"
        break
    counter +=1
print ""
print '{:>22s}'.format("number of iterations: "), counter
if converge ==1:
    rvector2 = EqtoEc(rvector2)
    rDotVector2 = EqtoEc(rDotVector2)
    r2 = np.linalg.norm(rvector2)
    a2, e2, i2, O2, w2, M2, E2 = babyOD(rvector2, rDotVector2)
    printrDotVector2 = k*rDotVector2
    print '{:>22s}'.format("range: "), rho1mag
    print '{:>22s}'.format("rvector2: "), rvector2
    print '{:>22s}'.format("rDotVector2: "), printrDotVector2
    print '{:>22s}'.format("a: "), a2
    print '{:>22s}'.format("e: "), e2
    print '{:>22s}'.format("i: "), degrees(i2)
    print '{:>22s}'.format("O: "), degrees(O2)
    print '{:>22s}'.format("w: "), degrees(w2)
    print '{:>22s}'.format("M: "), degrees(M2)

################################## Differential Corrections #######################################
dc = input("\nCompute differential corrections?\n(1) Yes\n(0) No\nEnter number: ")
if dc == 1:
    nyears= []
    nmonths = []
    ndays = []
    nhms = []
    nRAs = []
    nDecs = []
    nRvectors = []

    for item in inputFile:
        item = list(item)
        nyears.append(item[0])
        nmonths.append(item[1])
        ndays.append(item[2])
        nhms.append([int(number) for number in item[3].split(':')])
        nRAs.append(item[4:7])
        nDecs.append(item[7:10])
        nRvectors.append(ast.literal_eval(item[10]))

    nhr1 = nhms[0][0]
    nhr2 = nhms[1][0]
    nhr3 = nhms[2][0]
    nhr4 = nhms[3][0]

    nMin1 = nhms[0][1]
    nMin2 = nhms[1][1]
    nMin3 = nhms[2][1]
    nMin4 = nhms[3][1]

    nsec1 = nhms[0][2]
    nsec2 = nhms[1][2]
    nsec3 = nhms[2][2]
    nsec4 = nhms[3][2]

    nmon1 = nmonths[0]
    nmon2 = nmonths[1]
    nmon3 = nmonths[2]
    nmon4 = nmonths[3]

    nday1 = ndays[0]
    nday2 = ndays[1]
    nday3 = ndays[2]
    nday4 = ndays[3]

    nyr1 = nyears[0]
    nyr2 = nyears[1]
    nyr3 = nyears[2]
    nyr4 = nyears[3]

    nRAhr1 = nRAs[0][0]
    nRAhr2 = nRAs[1][0]
    nRAhr3 = nRAs[2][0]
    nRAhr4 = nRAs[3][0]

    nRAMin1 = nRAs[0][1]
    nRAMin2 = nRAs[1][1]
    nRAMin3 = nRAs[2][1]
    nRAMin4 = nRAs[3][1]

    nRAsec1 = nRAs[0][2]
    nRAsec2 = nRAs[1][2]
    nRAsec3 = nRAs[2][2]
    nRAsec4 = nRAs[3][2]

    ndeg1 = nDecs[0][0]
    ndeg2 = nDecs[1][0]
    ndeg3 = nDecs[2][0]
    ndeg4 = nDecs[3][0]
    ##print ndeg1, ndeg2, ndeg3, ndeg4 

    namin1 = nDecs[0][1]
    namin2 = nDecs[1][1]
    namin3 = nDecs[2][1]
    namin4 = nDecs[3][1]
    ##print nmin1, ndeg2, ndeg3, ndeg4 

    nasec1 = nDecs[0][2]
    nasec2 = nDecs[1][2]
    nasec3 = nDecs[2][2]
    nasec4 = nDecs[3][2]
    ##print ndeg1, ndeg2, ndeg3, ndeg4 

    nRvector1 = nRvectors[0]
    nRvector2 = nRvectors[1]
    nRvector3 = nRvectors[2]
    nRvector4 = nRvectors[3]
    nRvector = [nRvector1, nRvector2, nRvector3, nRvector4]

    nt1 = CDtoJD(nhr1, nMin1, nsec1, nyr1, nmon1, nday1)
    nt2 = CDtoJD(nhr2, nMin2, nsec2, nyr2, nmon2, nday2)
    nt3 = CDtoJD(nhr3, nMin3, nsec3, nyr3, nmon3, nday3)
    nt4 = CDtoJD(nhr4, nMin4, nsec4, nyr4, nmon4, nday4)
    time = [nt1, nt2, nt3, nt4]
    timemid = 2457956.75 #time for epoch july 22, 2017 6 UT

    nRA1 = radRA(nRAhr1, nRAMin1, nRAsec1)
    nRA2 = radRA(nRAhr2, nRAMin2, nRAsec2)
    nRA3 = radRA(nRAhr3, nRAMin3, nRAsec3)
    nRA4 = radRA(nRAhr4, nRAMin4, nRAsec4)
    nRA = [nRA1, nRA2, nRA3, nRA4]

    nDec1 = radDec(ndeg1, namin1, nasec1)
    nDec2 = radDec(ndeg2, namin2, nasec2)
    nDec3 = radDec(ndeg3, namin3, nasec3)
    nDec4 = radDec(ndeg4, namin4, nasec4)
##    print nDec1, nDec2, nDec3, nDec4 


    calcRA1, calcDec1 = ephem(a2, e2, i2, O2, w2, M2, time[0], timemid, nRvector[0])
    calcRA2, calcDec2 = ephem(a2, e2, i2, O2, w2, M2, time[1], timemid, nRvector[1])
    calcRA3, calcDec3 = ephem(a2, e2, i2, O2, w2, M2, time[2], timemid, nRvector[2])
    calcRA4, calcDec4 = ephem(a2, e2, i2, O2, w2, M2, time[3], timemid, nRvector[3])
##    print calcDec1, calcDec2, calcDec3, calcDec4 

    diffRA1 = nRA1 - calcRA1
    diffRA2 = nRA2 - calcRA2
    diffRA3 = nRA3 - calcRA3
    diffRA4 = nRA4 - calcRA4
    dR = np.array([diffRA1, diffRA2, diffRA3, diffRA4])
    print "\ndiffRA: ", diffRA1, diffRA2, diffRA3, diffRA4

    diffDec1 = nDec1 - calcDec1
    diffDec2 = nDec2 - calcDec2
    diffDec3 = nDec3 - calcDec3
    diffDec4 = nDec4 - calcDec4
    dD = np.array([diffDec1, diffDec2, diffDec3, diffDec4])
    diffList = np.array([diffRA1,diffDec1, diffRA2,diffDec2, diffRA3,diffDec3, diffRA4, diffDec4])
    print "diffDec: ", diffDec1, diffDec2, diffDec3, diffDec4
    oldRMS = sqrt((diffRA1**2+diffRA2**2+diffRA3**2+diffRA4**2+diffDec1**2+diffDec2**2+diffDec3**2+diffDec4**2)/2)
    print ""
    print '{:>22s}'.format("old RMS: "), oldRMS

    deltavpos = 1+10**-4
    deltavneg = 1-10**-4
    deltav = 10**-4


    ddxp = []
    ddxn = []
    ddyp = []
    ddyn = []
    ddzp = []
    ddzn = []
    ddxdotp = []
    ddxdotn = []
    ddydotp = []
    ddydotn = []
    ddzdotp = []
    ddzdotn = []
    a2list = []
    for i in range(4):
        nrvector2 = np.copy(rvector2)
        nrvector2[0] = rvector2[0]+deltav
        a2, e2, i2, O2, w2, M2, E2 = babyOD(nrvector2, rDotVector2)
        newRA1, newDec1 = ephem(a2, e2, i2, O2, w2, M2, time[i], timemid, nRvector[i])
        ddxp.append(newRA1)
        ddxp.append(newDec1)
        
        nrvector2 = np.copy(rvector2)
        nrvector2[0] = rvector2[0]-deltav
        a2, e2, i2, O2, w2, M2, E2 = babyOD(nrvector2, rDotVector2)
        newRA1, newDec1 = ephem(a2, e2, i2, O2, w2, M2, time[i], timemid, nRvector[i])
        ddxn.append(newRA1)
        ddxn.append(newDec1)
        
        nrvector2 = np.copy(rvector2)
        nrvector2[1] = rvector2[1]+deltav
        a2, e2, i2, O2, w2, M2, E2 = babyOD(nrvector2, rDotVector2)
        newRA1, newDec1 = ephem(a2, e2, i2, O2, w2, M2, time[i], timemid, nRvector[i])
        ddyp.append(newRA1)
        ddyp.append(newDec1)
        
        nrvector2 = np.copy(rvector2)
        nrvector2[1] = rvector2[1]-deltav
        a2, e2, i2, O2, w2, M2, E2 = babyOD(nrvector2, rDotVector2)
        newRA1, newDec1 = ephem(a2, e2, i2, O2, w2, M2, time[i], timemid, nRvector[i])
        ddyn.append(newRA1)
        ddyn.append(newDec1)

        nrvector2 = np.copy(rvector2)
        nrvector2[2] = rvector2[2]+deltav
        a2, e2, i2, O2, w2, M2, E2 = babyOD(nrvector2, rDotVector2)
        newRA1, newDec1 = ephem(a2, e2, i2, O2, w2, M2, time[i], timemid, nRvector[i])
        ddzp.append(newRA1)
        ddzp.append(newDec1)

        nrvector2 = np.copy(rvector2)
        nrvector2[2] = rvector2[2]-deltav
        a2, e2, i2, O2, w2, M2, E2 = babyOD(nrvector2, rDotVector2)
        newRA1, newDec1 = ephem(a2, e2, i2, O2, w2, M2, time[i], timemid, nRvector[i])
        ddzn.append(newRA1)
        ddzn.append(newDec1)


##        nrvector2 = np.copy(rvector2) # reset nrvector2

        
        nrDotVector2 = np.copy(rDotVector2)
        nrDotVector2[0] = nrDotVector2[0]+deltav
        a2, e2, i2, O2, w2, M2, E2 = babyOD(rvector2, nrDotVector2)
        newRA1, newDec1 = ephem(a2, e2, i2, O2, w2, M2, time[i], timemid, nRvector[i])
        ddxdotp.append(newRA1)
        ddxdotp.append(newDec1)

        nrDotVector2 = np.copy(rDotVector2)
        nrDotVector2[0] = nrDotVector2[0]-deltav
        a2, e2, i2, O2, w2, M2, E2 = babyOD(rvector2, nrDotVector2)
        newRA1, newDec1 = ephem(a2, e2, i2, O2, w2, M2, time[i], timemid, nRvector[i])
        ddxdotn.append(newRA1)
        ddxdotn.append(newDec1)

        nrvector2 = np.copy(rvector2)
        nrDotVector2 = np.copy(rDotVector2)
        nrDotVector2[1] = nrDotVector2[1]+deltav
        a2, e2, i2, O2, w2, M2, E2 = babyOD(rvector2, nrDotVector2)
        newRA1, newDec1 = ephem(a2, e2, i2, O2, w2, M2, time[i], timemid, nRvector[i])
        ddydotp.append(newRA1)
        ddydotp.append(newDec1)

        nrDotVector2 = np.copy(rDotVector2)
        nrDotVector2[1] = nrDotVector2[1]-deltav
        a2, e2, i2, O2, w2, M2, E2 = babyOD(rvector2, nrDotVector2)
        newRA1, newDec1 = ephem(a2, e2, i2, O2, w2, M2, time[i], timemid, nRvector[i])
        ddydotn.append(newRA1)
        ddydotn.append(newDec1)

        nrvector2 = np.copy(rvector2)
        nrDotVector2 = np.copy(rDotVector2)
        nrDotVector2[2] = nrDotVector2[2]+deltav
        a2, e2, i2, O2, w2, M2, E2 = babyOD(rvector2, nrDotVector2)
        newRA1, newDec1 = ephem(a2, e2, i2, O2, w2, M2, time[i], timemid, nRvector[i])
        ddzdotp.append(newRA1)
        ddzdotp.append(newDec1)

        nrDotVector2 = np.copy(rDotVector2)
        nrDotVector2[2] = nrDotVector2[2]-deltav
        a2, e2, i2, O2, w2, M2, E2 = babyOD(rvector2, nrDotVector2)
        newRA1, newDec1 = ephem(a2, e2, i2, O2, w2, M2, time[i], timemid, nRvector[i])
        ddzdotn.append(newRA1)
        ddzdotn.append(newDec1)
        a2list.append(a2)
##    print a2list
##    print a2
    ddxp = np.array(ddxp)   
    ddxn = np.array(ddxn)
    ddyp = np.array(ddyp)
    ddyn = np.array(ddyn)
    ddzp = np.array(ddzp)
    ddzn = np.array(ddzn)
    ddxdotp = np.array(ddxdotp)
    ddxdotn = np.array(ddxdotn)
    ddydotp = np.array(ddydotp)
    ddydotn = np.array(ddydotn)
    ddzdotp = np.array(ddzdotp)
    ddzdotn = np.array(ddzdotn)
    ddx = (ddxp-ddxn)/(2*deltav)
    ddy = (ddyp-ddyn)/(2*deltav)
    ddz = (ddzp-ddzn)/(2*deltav)
    ddxdot = (ddxdotp-ddxdotn)/(2*deltav)
    ddydot = (ddydotp-ddydotn)/(2*deltav)
    ddzdot = (ddzdotp-ddzdotn)/(2*deltav)

    ##print ddx
    ##print ddy
    ##print ddz
    ##print ddxdot
    ##print ddydot
    ##print ddzdot

    #Vertical column 
    m1 = np.array([ddx, ddy, ddz, ddxdot, ddydot, ddzdot])
##    print m1

##    vermatrix = np.array([[m1[0,0]*dR[0] + m1[0,2]*dR[1] + m1[0,4]*dR[2] + m1[0,6]*dR[3] + m1[0,1]*dD[0] + m1[0,3]*dD[1] + m1[0,5]*dD[2] + m1[0,7]*dD[3]],
##                         [m1[1,0]*dR[0] + m1[1,2]*dR[1] + m1[1,4]*dR[2] + m1[1,6]*dR[3] + m1[1,1]*dD[0] + m1[1,3]*dD[1] + m1[1,5]*dD[2] + m1[1,7]*dD[3]],
##                         [m1[2,0]*dR[0] + m1[2,2]*dR[1] + m1[2,4]*dR[2] + m1[2,6]*dR[3] + m1[2,1]*dD[0] + m1[2,3]*dD[1] + m1[2,5]*dD[2] + m1[2,7]*dD[3]],
##                         [m1[3,0]*dR[0] + m1[3,2]*dR[1] + m1[3,4]*dR[2] + m1[3,6]*dR[3] + m1[3,1]*dD[0] + m1[3,3]*dD[1] + m1[3,5]*dD[2] + m1[3,7]*dD[3]],
##                         [m1[4,0]*dR[0] + m1[4,2]*dR[1] + m1[4,4]*dR[2] + m1[4,6]*dR[3] + m1[4,1]*dD[0] + m1[4,3]*dD[1] + m1[4,5]*dD[2] + m1[4,7]*dD[3]],
##                         [m1[5,0]*dR[0] + m1[5,2]*dR[1] + m1[5,4]*dR[2] + m1[5,6]*dR[3] + m1[5,1]*dD[0] + m1[5,3]*dD[1] + m1[5,5]*dD[2] + m1[5,7]*dD[3]]])
    vermatrix = np.array([[np.dot(ddx, diffList)],
                          [np.dot(ddy, diffList)],
                          [np.dot(ddz, diffList)],
                          [np.dot(ddxdot, diffList)],
                          [np.dot(ddydot, diffList)],
                          [np.dot(ddzdot, diffList)]])
    ##print vermatrix
    
    Jacobmatrix = np.dot(m1, m1.T)
##    print Jacobmatrix

    solnmatrix = np.dot(np.linalg.inv(Jacobmatrix), vermatrix)
##    print solnmatrix

    newrvector2 = np.array([(solnmatrix[0]+rvector2[0])[0], (solnmatrix[1]+rvector2[1])[0], (solnmatrix[2]+rvector2[2])[0]])
    newrDotVector2 = np.array([(solnmatrix[3]+rDotVector2[0])[0], (solnmatrix[4]+rDotVector2[1])[0], (solnmatrix[5]+rDotVector2[2])[0]])
    newprintrDotVector2 = newrDotVector2*k
##    print newrvector2, newrDotVector2
    a2, e2, i2, O2, w2, M2, E2 = babyOD(newrvector2, newrDotVector2)
    
    calcRA1, calcDec1 = ephem(a2, e2, i2, O2, w2, M2, time[0], timemid, nRvector[0])
    calcRA2, calcDec2 = ephem(a2, e2, i2, O2, w2, M2, time[1], timemid, nRvector[1])
    calcRA3, calcDec3 = ephem(a2, e2, i2, O2, w2, M2, time[2], timemid, nRvector[2])
    calcRA4, calcDec4 = ephem(a2, e2, i2, O2, w2, M2, time[3], timemid, nRvector[3])
    ##print calcDec1, calcDec2, calcDec3, calcDec4 

    diffRA1 = nRA1 - calcRA1
    diffRA2 = nRA2 - calcRA2
    diffRA3 = nRA3 - calcRA3
    diffRA4 = nRA4 - calcRA4
    print "\ndiffRA: ", diffRA1, diffRA2, diffRA3, diffRA4

    diffDec1 = nDec1 - calcDec1
    diffDec2 = nDec2 - calcDec2
    diffDec3 = nDec3 - calcDec3
    diffDec4 = nDec4 - calcDec4
    print "diffDec: ", diffDec1, diffDec2, diffDec3, diffDec4

    RMS = sqrt((diffRA1**2+diffRA2**2+diffRA3**2+diffRA4**2+diffDec1**2+diffDec2**2+diffDec3**2+diffDec4**2)/2)
    print '{:>22s}'.format("new rvector2: "), newrvector2
    print '{:>22s}'.format("new rDotVector2: "), newprintrDotVector2
    print '{:>22s}'.format("new a: "), a2
    print '{:>22s}'.format("new e: "), e2
    print '{:>22s}'.format("new i: "), degrees(i2)
    print '{:>22s}'.format("new O: "), degrees(O2)
    print '{:>22s}'.format("new w: "), degrees(w2)
    print '{:>22s}'.format("new M: "), degrees(M2)
    print '{:>22s}'.format("new RMS: "), RMS
