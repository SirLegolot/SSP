from __future__ import division
import math

def magnitude(vec):
    squareSum = 0
    for num in vec:
        squareSum = squareSum + num**2
    mag = math.sqrt(squareSum)
    return mag

print(magnitude([1,1,1,1]))

def dot(vec1, vec2):
    dotProduct = 0
    for num in range(len(vec1)):
        dotProduct = dotProduct + vec1[num]*vec2[num]
    return dotProduct

print(dot([2,5,6],[3,7,8]))

def cross(vec1,vec2):
    elem1 = vec1[1]*vec2[2] - vec2[1]*vec1[2]
    elem2 = vec1[2]*vec2[0] - vec2[2]*vec1[0]
    elem3 = vec1[0]*vec2[1] - vec2[0]*vec1[1]
    crossProduct = [elem1, elem2, elem3]
    return crossProduct

print cross([1,0,0],[0,0,1])
