# starter code for exercise 0 on programming homework 2

import numpy as np
import copy

fruits = np.array([["Apple","Banana","Blueberry","Cherry"],
["Coconut","Grapefruit","Kumquat","Mango"],
["Nectarine","Orange","Tangerine","Pomegranate"],
["Lemon","Raspberry","Strawberry","Tomato"]])

#a
print "a"
print fruits[3,3]

#b
print "b"
print fruits[1:3,1:3]

#c
print "c"
print fruits[::2]

#d
print "d"
print fruits[-2:-4:-1,-2:-4:-1]

#e
print "e"
fruitscopy  = copy.copy(fruits)
fruitscopy[:,[0,3]] = fruits[:,[3,0]]
print fruitscopy

#f
print "f"
fruits[:4,:4] = "SLICED!"
print fruits
