# Warm Up
# Exercises
# 6/22/2017
# Neel Gandhi

from __future__ import division

def sumList(nums):
    sum = 0
    for num in nums:
        sum = sum + num
    return sum 

print "testing sumList"
print sumList([])              # expected output: 0
print sumList([3])             # expected output: 3
print sumList([1., 4.5, -3.2]) # expected output: 2.3

def estimatePi(n):
    piapprox = 0
    for term in range(n):
        piapprox = piapprox + ((-1)**term) * 4/(2*term +1)
    return piapprox

print "testing estimatePi"
print estimatePi(1)     # expected (approximate) output: 4.0
print estimatePi(10)    # expected (approximate) output: 3.0418396189294032
print estimatePi(100)   # expected (approximate) output: 3.1315929035585537
print estimatePi(1000)  # expected (approximate) output: 3.140592653839794
print estimatePi(10000) # expected (approximate) output: 3.1414926535900345

def scaleVec(vec, scalar):
    scaledList = []
    for num in vec:
        scaledList.append(num*scalar)
    print scaledList

print "testing scaleVec"
vec = []
print scaleVec(vec, 5) # expected output: []
print vec # expected output: []
vec = [1]
scaleVec(vec, 5) # expected output: [5]
print vec # expected output: [1]
vec = [-2, 1.5, 0.]
scaleVec(vec, 6) # expected output: [-12., 3., 0.]
print vec # expected output: [-2, 1.5, 0.]
