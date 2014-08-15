"""
This program plots the 'Soil Heat Flux' for two soil types in
Athens, Ga. The surface temperature is read from a .csv file and is
taken as the monthly average of ambient air temperature. Surface
temperature at 3 m below the surface is taken as the yearly average
of ambient air temperatures.
"""

import numpy as np
from numpy import genfromtxt
import matplotlib.pyplot as plt

z = 3 # m
tz = 0 # K
kls = 0.113 # W / m / K
kws = 2.68 # W / m / K

temps = np.zeros(12)
lsoil = np.zeros(12)
ssoil = np.zeros(12)
x = np.zeros(12)

input_data = genfromtxt('athensTemp.csv', delimiter=',')
data = input_data[1:13,1]

for i in range(len(data)):
    temps[i] = data[i]
    x[i] = i + 1
    tz += data[i]

avg_tz = tz / 12

for i in range(len(data)):
    lsoil[i] = kls * (temps[i] - avg_tz) / z
    ssoil[i] = kws * (temps[i]- avg_tz) / z

#######################################################
### PLOT CONFIG 

plt.figure(1)
#plt.subplot(211)
#plt.axis([0,13,-17,3])
plt.title('Rate of Heat Transfer from or into the Soil')
plt.scatter(x,lsoil, label = 'light w/ roots', s=35)
plt.legend(loc=8)
#plt.subplot(212)
plt.scatter(x, ssoil, color='blue', marker = '+', s = 100, label = 'wet, sandy')
plt.legend(loc=8)
plt.xlabel('Month of Year')
plt.ylabel('Soil Heat Flux (W/m^2)')
plt.axis([0,13,-12,12])
plt.show()
### END OF PLOT CONFIG

#######################################################
### DISCUSSION
"""
The greatest heat flux happens during the hottest times of the year--summer. The soil that is more saturated with water shows greater variation in heat flux. 
"""
#######################################################
