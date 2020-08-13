# Exercies 1
# Exercises
# 6/22/2017
# Neel Gandhi

from __future__ import division
import math
from math import radians, pi

degrees = input("Enter degrees: ")
minutes = input("Enter minutes: ")
seconds = input("Enter seconds: ")
convToRad = raw_input("Convert to radians? ")
normalize = raw_input("Normalize? ")

if degrees<0:
    minutes = -minutes
    seconds = -seconds

def convert(deg, minu, sec):
    total = deg + minu/60 + sec/3600
    if convToRad == "y":
        total = total * math.pi /180
    if normalize == "y":
        while total<0:
            total = total +360
        while total>360:
            total = total -360
    return total

print(convert(degrees, minutes, seconds))

# GET ADDITIONAL DATA FROM USER HERE
# hint: use input function for numerical input and raw_input function for text input

# perform angle conversion
# YOUR CODE HERE

# print results
# YOUR CODE HERE

# test cases ('y' and 'n' indicate response to either converting to radians or to normalization, converting to radians given first)
# if the user inputs 90 6 36 n n, your program should print 90.11
# if the user inputs 90 6 36 y y, your program should print 1.57271618897
# if the user inputs -90 6 36 y n, your program should print -1.57271618897
# if the user inputs -90 6 36 n y, your program should print 269.89
# if the user inputs 540 0 0 n y, your program should print 180.0
