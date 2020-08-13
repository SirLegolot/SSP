#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Orbit Visualization Program
# written by Michael Kelessoglou SSP '10

# The purpose of this program is to acquaint the untrained user with the orbital elements. A target Near Earth
# Asteroid is placed in a model Solar System. The user has the ability to modify the asteroid's orbital elements,
# as well as the speed of the animation. The user is provided with instructions of how to use the controls and
# information about the orbital elements



# built-in functions, objects, and other properties are loaded from modules

from visual import *
from visual.controls import *
from math import *



# function converts meters to astronomical units

def AU(m):
    au = m/(149.60e9)
    return au


# function rotates a celestial body and its orbit representation (and its orbital elements representations if it is the asteroid)
# from the ecliptic to its orbital plane

def rotateframe(body):
    body.f = frame()
    planevectory = arrow(opacity=0, pos=(0,0,0), axis=(0,1,0), frame=body.f)
    planevectorz = arrow(opacity=0, pos=(0,0,0), axis=(0,0,-1), frame=body.f)
    body.f.rotate(angle=body.Omega-3*pi/2, axis=planevectorz.axis, origin = (0,0,0))
    body.f.rotate(angle=body.i, axis=planevectory.axis, origin = (0,0,0))
    body.f.rotate(angle=body.w-pi/2, axis=planevectorz.axis, origin = (0,0,0))
    if body == asteroid:
        semimajoraxis.pos = (-body.a*body.e,0,0)
        semimajoraxis.axis = (body.a,0,0)
        ascendingnode.axis = (0,-body.a-1,0)
        periheliumradius.axis = (body.a-body.a*body.e,0,0)
        i1.axis = (body.a+1,0,0)
        planevectory = arrow(opacity=0, pos=(0,0,0), axis=(0,-1,0))
        planevectorz = arrow(opacity=0, pos=(0,0,0), axis=(0,0,1))
        planevectory.rotate(angle=body.Omega-3*pi/2, axis=planevectorz.axis, origin = (0,0,0))
        semimajoraxis.rotate(angle=body.Omega-3*pi/2, axis=planevectorz.axis, origin = (0,0,0))
        ascendingnode.rotate(angle=body.Omega-3*pi/2, axis=planevectorz.axis, origin = (0,0,0))
        periheliumradius.rotate(angle=body.Omega-3*pi/2, axis=planevectorz.axis, origin = (0,0,0))
        i1.rotate(angle=body.Omega-3*pi/2, axis=planevectorz.axis, origin = (0,0,0))
        planevectorz.rotate(angle=body.i, axis=planevectory.axis, origin = (0,0,0))
        semimajoraxis.rotate(angle=body.i, axis=planevectory.axis, origin = (0,0,0))
        periheliumradius.rotate(angle=body.i, axis=planevectory.axis, origin = (0,0,0))
        i1.rotate(angle=body.i, axis=planevectory.axis, origin = (0,0,0))
        semimajoraxis.rotate(angle=body.w-pi/2, axis=planevectorz.axis, origin = (0,0,0))
        periheliumradius.rotate(angle=body.w-pi/2, axis=planevectorz.axis, origin = (0,0,0))
        focallength.pos = semimajoraxis.pos
        focallength.axis = semimajoraxis.axis*body.e
        i2.axis = (i1.axis.x,i1.axis.y,0)


# function creates a display of the celestial body's orbit
        
def orbit(body):
    ellipse = curve(color=body.color, radius=0.002, frame = body.frame, display=SolarSystem)
    theta = 0.
    while theta < 2*pi+0.001:
        world_pos = (body.a*cos(theta)-body.a*body.e,sqrt(body.a**2-(body.a*body.e)**2)*sin(theta),0.)
        ellipse.append(body.f.world_to_frame(world_pos))
        theta = theta + 0.001
    return ellipse


# function tracks the celestial body on its orbit

def track(body):
    E = 0.
    #Newton-Raphson method inverts the equation and gives very accurate solution
    count = 0
    while abs(E - body.e*sin(E) - body.M) > 0.000001:
        E = E - (E - body.e*sin(E) - body.M)/(1 - body.e*cos(E))
        count = count + 1
        if count > 1000:
            break #stops the loop if very accurate solution cannot be found, so that animation does not stop
    if body == asteroid:
        if abs(globals()["Ememory"] - E) > 1.: #repeats Newton-Raphson method with different initial value if first attempt did not converge
            E = pi
            count = 0
            while abs(E - body.e*sin(E) - body.M) > 0.001:
                E = E - (E - body.e*sin(E) - body.M)/(1 - body.e*cos(E))
                count = count + 1
                if count > 100:
                    break
        globals()["Ememory"] = E #helps detect future errors in the Newton-Raphson method
    world_pos = (body.a*(cos(E) - body.e),sqrt(body.a**2 - (body.a*body.e)**2)*sin(E),0.)
    body.pos = body.f.world_to_frame(world_pos)
    body.label.pos = body.pos

        
# functions that operate the controls

def set_semimajor_axis(a):
    asteroid.semimajoraxis = a.value

def set_eccentricity(e):
    asteroid.eccentricity = e.value

def set_inclination(i):
    asteroid.inclination = i.value

def set_longitude_of_ascending_node(Omega):
    asteroid.longitude_of_ascending_node = Omega.value

def set_argument_of_perihelium(w):
    asteroid.argument_of_perihelium = w.value

def set_animation_speed(speed):
    asteroid.animation_speed = speed.value

def set_body_size():
    if size_toggle.value:
        count = 0
        while count < 8:
            planets[count].radius = planets[count].trueradius*1000
            count = count + 1
        asteroid.radius = asteroid.trueradius*500000
        sun.radius = sun.trueradius*10
    else:
        count = 0
        while count < 8:
            planets[count].radius = planets[count].trueradius
            count = count + 1
        asteroid.radius = asteroid.trueradius
        sun.radius = sun.trueradius

def set_label_visibility():
    if label_toggle.value:
        count = 0
        while count < 8:
            planets[count].label.box = False
            planets[count].label.line = False
            planets[count].label.text = ''
            planets[count].label.opacity = 0
            count = count + 1
        asteroid.label.box = False
        asteroid.label.line = False
        asteroid.label.text = ''
        asteroid.label.opacity = 0
    else:
        count = 0
        while count < 8:
            planets[count].label.box = True
            planets[count].label.line = True
            planets[count].label.opacity = 1
            count = count + 1
        mercury.label.text = 'Mercury'
        venus.label.text = 'Venus'
        earth.label.text = 'Earth'
        mars.label.text = 'Mars'
        jupiter.label.text = 'Jupiter'
        saturn.label.text = 'Saturn'
        uranus.label.text = 'Uranus'
        neptune.label.text = 'Neptune'
        asteroid.label.box = True
        asteroid.label.line = True
        asteroid.label.text = 'Asteroid'
        asteroid.label.opacity = 1

def continueinstructions(number):
    number = number + 1
    if number < 23:
        instructionstext.text = instructions[number]
        globals()["instructionsnumber"] = number

def backinstructions(number):
    if number > 0:
        number = number - 1
        instructionstext.text = instructions[number]
        globals()["instructionsnumber"] = number

def showa():
    semimajoraxis.opacity = 1
    focallength.opacity = 0
    vernalequinox.opacity = 0
    ascendingnode.opacity = 0
    periheliumradius.opacity = 0
    i1.opacity = 0
    i2.opacity = 0

def showe():
    semimajoraxis.opacity = 1
    focallength.opacity = 1
    vernalequinox.opacity = 0
    ascendingnode.opacity = 0
    periheliumradius.opacity = 0
    i1.opacity = 0
    i2.opacity = 0

def showi():
    semimajoraxis.opacity = 0
    focallength.opacity = 0
    vernalequinox.opacity = 0
    ascendingnode.opacity = 0
    periheliumradius.opacity = 0
    i1.opacity = 1
    i2.opacity = 1

def showOmega():
    semimajoraxis.opacity = 0
    focallength.opacity = 0
    vernalequinox.opacity = 1
    ascendingnode.opacity = 1
    periheliumradius.opacity = 0
    i1.opacity = 0
    i2.opacity = 0

def showw():
    semimajoraxis.opacity = 0
    focallength.opacity = 0
    vernalequinox.opacity = 0
    ascendingnode.opacity = 1
    periheliumradius.opacity = 1
    i1.opacity = 0
    i2.opacity = 0

def clear():
    semimajoraxis.opacity = 0
    focallength.opacity = 0
    vernalequinox.opacity = 0
    ascendingnode.opacity = 0
    periheliumradius.opacity = 0
    i1.opacity = 0
    i2.opacity = 0



# Display for the simulation
SolarSystem = display(x=0, y=0, width=700, height=700, range=5, forward=vector(0,0,-1), newzoom=1)


# Creation of the sun, the planets, and the asteroid (positions are accurate as of 01/01/00)

sun = sphere(pos=(0.,0.,0.), radius=AU(6.96e8), trueradius=AU(6.96e8), color=color.yellow, light=local_light(pos=(0.,0.,0.), color=color.white), material=materials.emissive)
mercury = sphere(a=0.387098, e=0.205630, i=radians(7.005), Omega=radians(48.331), w=radians(29.124), M=radians(174.796), radius=AU(2440.e3), trueradius=AU(2440.e3), color=(0.888,0.736,0.540))
venus = sphere(a=0.723332, e=0.0068, i=radians(3.39471), Omega=radians(76.67069), w=radians(54.85229), M=radians(50.44675), radius=AU(6051.8e3), trueradius=AU(6051.8e3), color=color.magenta)
earth = sphere(a=1.00000261, e=0.01671123, i=0., Omega=radians(348.73936), w=radians(114.20783), M=radians(357.51716), radius=AU(6371.01e3), trueradius=AU(6371.01e3), color=color.blue)
mars = sphere(a=1.523679, e=0.093315, i=radians(1.850), Omega=radians(49.562), w=radians(286.537), M=radians(19.3564), radius=AU(3389.9e3), trueradius=AU(3389.9e3), color=color.red)
jupiter = sphere(a=5.204267, e=0.048775, i=radians(1.305), Omega=radians(100.492), w=radians(275.066), M=radians(18.818), radius=AU(69911e3), trueradius=AU(69911e3), color=color.orange)
saturn = sphere(a=9.58201720, e=0.055723219, i=radians(2.485240), Omega=radians(113.642811), w=radians(336.013862), M=radians(320.346750), radius=AU(58232e3), trueradius=AU(58232e3), color=color.yellow)
uranus = sphere(a=19.22941195, e=0.044405586, i=radians(0.772556), Omega=radians(73.989821), w=radians(96.541318), M=radians(142.955717), radius=AU(25362e3), trueradius=AU(25362e3), color=color.cyan)
neptune = sphere(a=30.10366151, e=0.011214269, i=radians(1.767975), Omega=radians(131.794310), w=radians(265.646853), M=radians(267.767281), radius=AU(24624e3), trueradius=AU(24624e3), color=color.blue)
asteroid = sphere(a=2.014, e=0.3897, i=radians(1.282), Omega=radians(280.88), w=radians(9.095), M=radians(312.6623), radius=AU(10e3), trueradius=AU(10e3))#initial orbital elements are those of 2000 NF5
planets = [mercury, venus, earth, mars, jupiter, saturn, uranus, neptune]


# visual representations of orbital elements

semimajoraxis = cylinder(pos=(-asteroid.a*asteroid.e,0,0), axis=(asteroid.a,0,0), radius=0.006, opacity=0)
vernalequinox = cylinder(pos=(0,0,0), axis=(asteroid.a+1,0,0), radius=0.006, opacity=0)
ascendingnode = cylinder(pos=(0,0,0), axis=((asteroid.a+1)*cos(asteroid.Omega),(asteroid.a+1)*sin(asteroid.Omega),0), radius=0.006, opacity=0)
periheliumradius = cylinder(pos=(0,0,0), axis=(asteroid.a-asteroid.a*asteroid.e,0,0), radius=0.006, opacity=0)
i1 = cylinder(pos=(0,0,0), axis=(asteroid.a+1,0,0), radius=0.006, opacity=0)
planevectory = arrow(opacity=0, pos=(0,0,0), axis=(0,1,0))
planevectorz = arrow(opacity=0, pos=(0,0,0), axis=(0,0,1))
semimajoraxis.rotate(angle=asteroid.Omega-3*pi/2, axis=planevectorz.axis, origin = (0,0,0))
periheliumradius.rotate(angle=asteroid.Omega-3*pi/2, axis=planevectorz.axis, origin = (0,0,0))
i1.rotate(angle=asteroid.Omega-3*pi/2, axis=planevectorz.axis, origin = (0,0,0))
planevectory.rotate(angle=asteroid.Omega-3*pi/2, axis=planevectorz.axis, origin = (0,0,0))
semimajoraxis.rotate(angle=asteroid.i, axis=planevectory.axis, origin = (0,0,0))
periheliumradius.rotate(angle=asteroid.i, axis=planevectory.axis, origin = (0,0,0))
i1.rotate(angle=asteroid.i, axis=planevectory.axis, origin = (0,0,0))
planevectorz.rotate(angle=asteroid.i, axis=planevectory.axis, origin = (0,0,0))
semimajoraxis.rotate(angle=asteroid.w-pi/2, axis=planevectorz.axis, origin = (0,0,0))
periheliumradius.rotate(angle=asteroid.w-pi/2, axis=planevectorz.axis, origin = (0,0,0))
focallength = cylinder(pos=semimajoraxis.pos, axis=semimajoraxis.axis*asteroid.e, radius=0.0065, color=color.red, opacity=0)
i2 = cylinder(pos=(0,0,0), axis=(i1.axis.x,i1.axis.y,0), radius=0.006, opacity=0)


# celestial body labels

mercury.label=label(pos=mercury.pos, xoffset=20, yoffset=10, border=2, height=15, text='Mercury')
venus.label=label(pos=venus.pos, xoffset=20, yoffset=10, border=2, height=15, text='Venus')
earth.label=label(pos=earth.pos, xoffset=20, yoffset=10, border=2, height=15, text='Earth')
mars.label=label(pos=mars.pos, xoffset=20, yoffset=10, border=2, height=15, text='Mars')
jupiter.label=label(pos=jupiter.pos, xoffset=20, yoffset=10, border=2, height=15, text='Jupiter')
saturn.label=label(pos=saturn.pos, xoffset=20, yoffset=10, border=2, height=15, text='Saturn')
uranus.label=label(pos=uranus.pos, xoffset=20, yoffset=10, border=2, height=15, text='Uranus')
neptune.label=label(pos=neptune.pos, xoffset=20, yoffset=10, border=2, height=15, text='Neptune')
asteroid.label=label(pos=asteroid.pos, xoffset=20, yoffset=10, border=2, height=15, text='Asteroid')


# tracks latest eccentric anomaly of asteroid

Ememory = 0. #its initial value might need to be changed if different defaults are given for the asteroid


# planets and asteroid are located and visual representations of orbits are created

count = 0
while count < 8:
    planets[count].f = frame()
    planets[count].frame = planets[count].f
    rotateframe(planets[count])
    planets[count].orbit = orbit(planets[count])
    track(planets[count])
    count = count + 1
asteroid.f = frame()
asteroid.frame = asteroid.f
rotateframe(asteroid)
asteroid.orbit = orbit(asteroid)
track(asteroid)



#Instructions window

instructionsdisplay = display(x=700, y=60, width=400, height=400, range=10)

instructions = [u"Welcome to SSP's orbital\nastronomy tutorial. The purpose of\nthis program is to acquaint you\nwith some of the concepts to be\ncovered in SSP. Click 'NEXT'\nto continue (it's right above).",
                u"In the window to the left, you\ncan see a scaled model of\nthe Solar System. The coloured\ncircles are graphic representations\nof orbits. The positions of the\nplanets are indicated by labels.\n(Click 'NEXT')",
                u"You probably have noticed at\nthis point that neither the\nplanets nor the Sun are visible,\nwhich shouldn't alarm you since the\nSun's radius is about 6500 times\nsmaller than that of Neptune's\norbit. (Click 'NEXT')",
                u"Select the display of the simulation.\nYou can zoom by holding down\nboth the left and right buttons on your\nmouse and moving the mouse, or\n(on a Mac) command-clicking and\ndragging. You can rotate by holding\ndown the right button, or (on a\nMac) option-clicking and dragging.\nThe Sun will always be at the\ncenter of the display. (Click 'NEXT')",
                u"Zoom in a little to the point\nwhere you can no longer see\nJupiter's orbit. Notice that\nbesides the planets an asteroid is\nalso labeled. We will be discussing\nthe orbit of this asteroid.\n(Click 'NEXT')",
                u"Objects in the Solar System\norbit around the Sun because of\nits gravitational pull. Their\norbits are elliptical with the Sun\nat one focus given that they\nare far enough from everything\nelse to avoid significant pulls.\n(Click 'NEXT')",
                u"An orbit is defined by five numbers\ncalled orbital elements. A sixth\nelement defines the position of the\nasteroid on the orbit at a given point\nin time. Below this text you can see\nthe values of the asteroid's orbital\nelements (excluding the time one).\n(Click 'NEXT')",
                u"The first orbital element is the\nsemi-major axis (a), which is half\nthe sum of the greatest and the\nsmallest distances from the Sun to\nthe object in the orbit. It is usually\nmeasured in astronomical units, the\naverage Earth to Sun distance.\n(Click 'NEXT')",
                u"Make sure you have a good view of\nthe asteroid's orbit and change the\nsemi-major axis, using the slider\nnext to the 'a' button. Watch how\nthe orbit changes. You can press the\n'a' button to create a representation\nof the semi-major axis to help you.\n(Click 'NEXT')",
                u"If labels are getting in the way of your\nview, you can turn them off by setting the\ncorresponding switch to 'Labels OFF'.\nIf you would like to see the planets\nand the asteroid, you can compromise\nscaling accuracy by setting the other\nswitch to 'Visible'. Click on the lever\nto change a switch. (Click 'NEXT')",
                u"The second orbital element is\neccentricity, the ratio of the\nfocal length of the orbit's ellipse to\nits semi-major axis. The focal\nlength is the distance from the\ncenter of the ellipse to the focus\n(Sun). Eccentricity varies between\n0 and 1. (Click 'NEXT')",
                u"Change the eccentricity of the orbit,\nusing the corresponding slider. The\nbutton will create visual\nrepresentations of the semi-major\naxis (white) and the focal length\n(red), to help you understand what\nyou are changing.\n(Click 'NEXT')",
                u"The first two orbital elements\ndefine the size and shape of\nthe orbit, but say nothing about its\norientation in space. Make sure you\nrotate so that you don't have an\noverhead view of the orbit.\n(Click 'NEXT')",
                u"The third orbital element is\ninclination(i), the angle between the\nplane of the orbit and the ecliptic,\nthe plane of the Earth's orbit. Alter\nthe inclination using its slider. The\nbutton creates two lines, the angle\nbetween which is the inclination.\n(Click 'NEXT')",
                u"The vernal point is the direction in\nthe sky (celestial sphere), through\nwhich the Sun passes on the spring\nequinox on its way from south to\nnorth. The ascending node is the\npoint in the orbit at which the\nobject rises above the ecliptic.\n(Click 'NEXT')",
                u"The fourth orbital element is the\nlongitude of the ascending node(Î©),\nthe angle between the vernal point\nand the ascending node. Change it\nusing its slider and observe. The\nbutton creates representations of the\nvernal point and the ascending node.\n(Click 'NEXT')",
                u"The fifth orbital element is the\nargument of the perihelium (Ï‰), the\nangle between the ascending node\nand the perihelium point, the point\nin the orbit at which the object is\nclosest to the Sun.\n(Click 'NEXT')",
                u"Change the argument of the\nperihelium using the corresponding\nslider and watch how the orbit\nchanges. Notice that the plane of the\norbit is already defined. The button\ncreates visual representations of\nthe ascending node and perihelium.\n(Click 'NEXT')",
                u"The sixth orbital element manifests\nitself in many forms. One is the\ntime of perihelium passage (T), a\ntime at which the object passed\nfrom the perihelium. Another is the\nmean anomaly(M), an angle that\nchanges linearly with time.\n(Click 'NEXT')",
                u"The bottom slider controls the\nanimation speed. Slide it and watch\nthe movement(you need to have\neither 'Visible' or 'Labels ON'\nselected). The mean anomaly\nconstantly changes, whereas T is\nconstant.\n(Click 'NEXT')",
                u"To clear the visual representations\nof orbital elements, click 'clear'.\nYou can play around a bit more if\nyou want, but we hope you\nlearned something about orbits\ntoday. To move to previous\ninstructions press 'PREVIOUS'.\n(Click 'NEXT')",
                u"Remember that the main project of\nSSP is the Orbital Determination,\nwhich uses telescope observations\nof an asteroid to determine its\norbital elements. The default orbital\nelements were the ones my team\ngot for our asteroid, 2000 NF5.\n(Click 'NEXT')",
                u"This tutorial is over.\nThank you for completing it and\nwe hope to see you at SSP."]
instructionsnumber = 0
instructionstext = label(pos=(0,4,0), xoffset=0, yoffset=0, box=False, line=False, height=22, text = instructions[instructionsnumber], opacity=0, font='serif')


# text displaying values of orbital elements of target asteroid

atext = label(pos=(-6,-3), xoffset=0, yoffset=0, box=False, line=False, height=24, opacity=0, font='serif')
etext = label(pos=(-7.2,-5), xoffset=0, yoffset=0, box=False, line=False, height=24, opacity=0, font='serif')
itext = label(pos=(-7,-7), xoffset=0, yoffset=0, box=False, line=False, height=24, opacity=0, font='serif')
Omegatext = label(pos=(4,-3), xoffset=0, yoffset=0, box=False, line=False, height=24, opacity=0, font='serif')
wtext = label(pos=(3.5,-4.8), xoffset=0, yoffset=0, box=False, line=False, height=24, opacity=0, font='serif')
speedtext = label(pos=(3,-7.2), xoffset=0, yoffset=0, box=False, line=False, height=24, opacity=0, font='serif')



# Controls display

ControlPanel = controls(x=700, y=460, width=400, height=240, range=200)

a_slider = slider(pos=(-140,100), width=10, length=300, axis=(1,0,0), min=0.00001, max=5., action=lambda: set_semimajor_axis(a_slider))
e_slider = slider(pos=(-140,70), width=10, length=300, axis=(1,0,0), min=0., max=0.999, action=lambda: set_eccentricity(e_slider))
i_slider = slider(pos=(-140,40), width=10, length=300, axis=(1,0,0), min=0., max=pi/2, action=lambda: set_inclination(i_slider))
Omega_slider = slider(pos=(-140,10), width=10, length=300, axis=(1,0,0), min=0., max=2*pi, action=lambda: set_longitude_of_ascending_node(Omega_slider))
w_slider = slider(pos=(-140,-20), width=10, length=300, axis=(1,0,0), min=0., max=2*pi, action=lambda: set_argument_of_perihelium(w_slider))
speed_slider = slider(pos=(-140,-95), width=10, length=300, axis=(1,0,0), min=0., max=2000., action=lambda: set_animation_speed(speed_slider))
size_toggle = toggle(pos=(-50,-55), width=18, height=18, text0='Real Size', text1='Visible', action=lambda: set_body_size())
label_toggle = toggle(pos=(50,-55), width=18, height=18, text0='Labels ON', text1='Labels OFF', action=lambda: set_label_visibility())
abutton = button(pos=a_slider.pos-(20,0), width=20, height=20, text='a', action=lambda: showa())
ebutton = button(pos=e_slider.pos-(20,0), width=20, height=20, text='e', action=lambda: showe())
ibutton = button(pos=i_slider.pos-(20,0), width=20, height=20, text='i', action=lambda: showi())
Omegabutton = button(pos=Omega_slider.pos-(20,0), width=20, height=20, text=u'Î©', action=lambda: showOmega())
wbutton = button(pos=w_slider.pos-(20,0), width=20, height=20, text=u'Ï‰', action=lambda: showw())
clearbutton = button(pos=(-170,-50), width=50, height=20, text='clear', action=lambda: clear())


# initial values for sliders set

a_slider.value = asteroid.a
e_slider.value = asteroid.e
i_slider.value = asteroid.i
Omega_slider.value = asteroid.Omega
w_slider.value = asteroid.w
speed_slider.value = 0.

set_semimajor_axis(a_slider)
set_eccentricity(e_slider)
set_inclination(i_slider)
set_longitude_of_ascending_node(Omega_slider)
set_argument_of_perihelium(w_slider)
set_animation_speed(speed_slider)


# Second controls display

ControlPanel2 = controls(x=700, y=40, width=400, height=60, range=200)

instructionsbutton = button(pos=(70,0), width=90, height=40, text='NEXT', action=lambda: continueinstructions(instructionsnumber))
instructionsbackbutton = button(pos=(-70,0), width=90, height=40, text='PREVIOUS', action=lambda: backinstructions(instructionsnumber))



# initial displays of orbital elements set

atext.text = u"a = %s AU" % (round(asteroid.a,2))
etext.text = u"e = %s" % (round(asteroid.e,3))
itext.text = u"i = %sÂ°" % (round(degrees(asteroid.i),1))
Omegatext.text = u"Î© = %sÂ°" % (round(degrees(asteroid.Omega),1))
wtext.text = u"Ï‰ = %sÂ°" % (round(degrees(asteroid.w),1))
speedtext.text = "animation speed: \n%s million x real time" % (round(speed_slider.value*50000./1000000.,2))



# helper variables used during animation

check = 0
count2 = 0.
animationspeed = speed_slider.value


# The animation begins! (actually it begins after you change the speed_slider)

while True:
    rate(2000) #sets a maximum of 2000 loops per second to make sure animation will run at the same speed on all computers
    ControlPanel.interact()
    ControlPanel2.interact()
    # indicates text must be changed because the animation speed has been changed
    if not animationspeed == speed_slider.value:
        check = 1
        count2 = 0.
        animationspeed = speed_slider.value
    if speed_slider.value > 0.: #simulation can speed up to 100000000 times faster than real time
        count = 0
        while count < 8:
            planets[count].M = (planets[count].M + sqrt(1/planets[count].a**3)*0.009955106388286*speed_slider.value/2000.) % (2*pi) #mean anomaly changes
            track(planets[count]) #movement of planets
            count = count + 1
        asteroid.M = (asteroid.M + sqrt(1/asteroid.a**3)*0.009955106388286*speed_slider.value/2000.) % (2*pi) #mean anomaly changes
        track(asteroid) #movement of asteroid
    count2 = count2 + 1./2000. #keeps track of real time
    # sets the value of orbital elements to that of their corresponding sliders if one of them is changed
    if not (a_slider.value == asteroid.a and e_slider.value == asteroid.e and i_slider.value == asteroid.i and Omega_slider.value == asteroid.Omega and w_slider.value == asteroid.w):
        asteroid.orbit.pos = []
        asteroid.a = a_slider.value
        asteroid.e = e_slider.value
        asteroid.i = i_slider.value
        asteroid.Omega = Omega_slider.value
        asteroid.w = w_slider.value
        rotateframe(asteroid)
        asteroid.orbit = orbit(asteroid)
        track(asteroid)
        # indicates text must be changed because the orbital elements have been changed
        check = 1
        count2 = 0.
    # adjusts the displayed values of orbital elements and animation speed once a change is made
    if check == 1 and count2 > 0.1:
        atext.text = u"a = %s AU" % (round(asteroid.a,2))
        etext.text = u"e = %s" % (round(asteroid.e,3))
        itext.text = u"i = %sÂ°" % (round(degrees(asteroid.i),1))
        Omegatext.text = u"Î© = %sÂ°" % (round(degrees(asteroid.Omega),1))
        wtext.text = u"Ï‰ = %sÂ°" % (round(degrees(asteroid.w),1))
        speedtext.text = "animation speed: \n%s million x real time" % (round(speed_slider.value*50000./1000000.,2))
        check = 0
        

