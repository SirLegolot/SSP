# exercise
# exercise 2
# 6/23/2017
# Neel Gandhi

from math import pi, radians, degrees, sin, cos, asin, acos, atan2

# DEFINE FUNCTION CONVERTING HA AND DEC TO ALTITUDE AND AZIMUTH HERE
# note: to have a function return two values, just do: return value1, value2
lat = radians(34.0727)
def altazi(HADegrees, decDegrees):
    HA = radians(HADegrees)
    dec = radians(decDegrees)
    altitude = degrees(asin(sin(lat)*sin(dec)+cos(lat)*cos(dec)*cos(HA)))
    #azimuth = degrees(asin(-cos(dec)*sin(HA)/cos(radians(altitude))))
    sinazi = -cos(dec)*sin(HA)/cos(radians(altitude))
    cosazi = (sin(dec) - sin(radians(altitude))*sin(lat))/(cos(radians(altitude))*cos(lat))
    azimuth = degrees(atan2(sinazi, cosazi))
    if azimuth < 0:
        azimuth= azimuth+ 360
    print "Your altitude is:", altitude, "degrees"
    print "Your azimuth is:", azimuth, "degrees"
    return altitude, azimuth


print "testing HA/Dec to Alt/Az"
altazi(54, 36)
altazi(204, 75)


# DEFINE FUNCTION CONVERTING EQUATORIAL TO ECLIPTIC HERE
# note: to have a function return two values, just do: return value1, value2
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
    sinlon = yprime/cos(radians(latitude))
    coslon = xprime/cos(radians(latitude))
    lon = degrees(atan2(sinlon,coslon))
    if lon <0:
        lon = lon +360
    #if RADegrees>180:
    #    lon = 360 - degrees(acos(xprime/cos(radians(latitude))))
    #else:
    #    lon = degrees(acos(xprime/cos(radians(latitude))))
    print "Your latitude is:", latitude, "degrees"
    print "Your longitude is:", lon, "degrees"
    return latitude, lon


    

print "testing equatorial to ecliptic"
# INCLUDE TEST CASES FROM HOMEWORK
longlat(250,36)
longlat(100.11,80)
longlat(345,-40)

