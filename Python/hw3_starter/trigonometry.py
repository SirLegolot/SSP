# PURPOSE OF THIS CODE
# PROJECT
# NAME
# DATE
from __future__ import division
from math import sin, cos, asin, acos, radians, degrees, pi

# a function to determine the quadrant of an angle based on its sine and cosine (in radians)
# returns the angle in the correct quadrant (in radians)
def findQuadrant(sine, cosine):
    if cosine > 0 and sine > 0: #1
        return asin(sine)

    if cosine < 0 and sine > 0: #2
        return acos(cosine)

    if cosine < 0 and sine < 0: #3
        return pi - asin(sine)

    if cosine > 0 and sine < 0: #4
        return 2*pi + asin(sine)

# a function that given the values (in radians) of two sides and the included angle of a spheical triangle
# returns the values of the remaining side and two angles (in radians)
def SAS(a, B, c):
    # YOUR CODE HERE (part a)
    b = acos(cos(a)*cos(c)+sin(a)*sin(c)*cos(B)) # calculates B
    sinA = sin(B)/sin(b)*sin(a)
    cosA = (cos(a)-cos(b)*cos(c))/(sin(b)*sin(c))
    A = degrees(findQuadrant(sinA, cosA)) # uses the findquadrant function to calculate A
    sinC = sin(B)/sin(b)*sin(c)
    cosC = (cos(c)-cos(a)*cos(b))/(sin(a)*sin(b))
    C = degrees(findQuadrant(sinC,cosC)) # uses the findquadrant function to calculate C
    b = degrees(b)
    
    return b, A, C

# YOUR CODE HERE (part b)
##a = input("what is a: ")
##B = input("what is B: ")
##c = input("what is c: ")
a = 106
arad = radians(a)
B = 114
Brad = radians(B)
c = 42
crad = radians(c)
b, A, C = SAS(arad, Brad, crad)
print "b=", b
print "A=", A
print "C=", C

