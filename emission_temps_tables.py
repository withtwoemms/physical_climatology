import math 
import numpy
import pickle
from prettytable import PrettyTable

sigma = 5.67e-8 # W / m^2 / K^4
L = 3.9e26 # W

planets = ["Pluto", "Jupiter", "Earth", "Mercury"]

mean_distances = {0 : 5.9e12, 1 : 7.78e11, 2 : 1.5e11, 3 : 5.8e10}
albedos = {0 : .4, 1 : .34, 2 : .3, 3 : .058}

solar_constants = numpy.zeros([len(planets)])
emiss_temps = numpy.zeros([len(planets)])

for i in range(len(planets)):
    solar_constants[i] = L / 4. / math.pi / mean_distances[i]**2
    emiss_temps[i] = (solar_constants[i] * (1 - albedos[i]) / 4. / sigma) ** (1. / 4.)
    

x = PrettyTable()

x.add_column("Planet", planets)
x.add_column("Solar Constant (W/m^2)", solar_constants)
x.add_column("Emission Temperature (K)", emiss_temps)

print '\n'+'Emission Temperatures of Selected Planets'+'\n',x
print '* EMISSION TEMPERATURE represents the surface temperature'+'\n'+ \
	    'necessary to balance absorbed solar radiation. A lower' +'\n'+ \
		'emission temperature calculated for Jupiter implies that'+'\n'+ \
		'some other radiation source is present. This source could'+'\n'+ \
		'be radiation emitted from the acceleration of particles by'+'\n'+ \
		"Jupiter's massive gravitational and resultant magnetic field"+'\n'

"""
Now what if the SOLAR LUMINOSITY was 30% less...
"""

L70percent = L * .7

for i in range(len(planets)):
    solar_constants[i] = L70percent / 4. / math.pi / mean_distances[i]**2
    emiss_temps[i] = (solar_constants[i] * (1 - albedos[i]) / 4. / sigma) ** (1. / 4.)
    

y = PrettyTable()

y.add_column("Planet", planets)
y.add_column("Solar Constant (W/m^2)", solar_constants)
y.add_column("Emission Temperature (K)", emiss_temps)
	
print '\n'+'Emission Temperatures if Solar Luminosity is 30% Less'+'\n',y

"""
Let's say planetary albedo was affected by the sudden disappreaence of clouds
plus the planet is a graybody with a 96% emissivity
"""

thirdofalbedos = {}
for i in range(len(albedos)):
    thirdofalbedos[i] = albedos[i] / 3
	
# since len(albedos) == len(planets)...

    solar_constants[i] = L / 4. / math.pi / mean_distances[i]**2
    emiss_temps[i] = (solar_constants[i] * (1 - thirdofalbedos[i]) / 4. / sigma / 0.96) ** (1. / 4.)

z = PrettyTable()

z.add_column("Planet", planets)
z.add_column("Solar Constant (W/m^2)", solar_constants)
z.add_column("Emission Temperature (K)", emiss_temps)

print '\n'+'Emission Temperatures if Planetary Albedos Reduced to 1/3'+'\n',z 
