from math import pi, radians, degrees, sin, cos, asin, acos

epsilon = radians(23.45)
def longlat(RADegrees, decDegrees):
    RA = radians(RADegrees)
    declination = radians(decDegrees)
    x = cos(RA)*cos(declination)
    y = sin(RA)*cos(declination)
    z = sin(declination)
    xprime = x
    yprime = y*cos(epsilon) + z*sin(epsilon)
    zprime = y*-sin(epsilon) + z*cos(epsilon)
    latitude = degrees(asin(zprime))
    lon = 360 - degrees(acos(xprime/cos(radians(latitude))))
    print "Your latitude is:", latitude, "degrees"
    print "Your longitude is:", lon, "degrees"
    return latitude, lon


    

print "testing equatorial to ecliptic"
# INCLUDE TEST CASES FROM HOMEWORK
longlat(345,-40)
