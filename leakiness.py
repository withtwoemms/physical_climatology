"""
This program plots thepy "leakiness" parameter, beta, as
a function of surface temperature, T_s.
"""

import numpy as np
import matplotlib.pyplot as plt

Se = 1379 # W / m^2
Sm = 597 # W / m^2
ae = 0.3 
am = 0.16
sigma = 5.67e-8 

interval = 100

#######################################################
### FILLING ARRAYS

beta1 = np.zeros([interval,2])
x1 = np.linspace(0.,1.,interval)

for i in range(interval):
    beta1[i,0] = x1[i]
    beta1[i,1] = ((Se * (1 - ae)) / 4 / sigma / (1 - (1 - x1[i]) / 2)) ** (1./4.)

#------------------------------------------------------

beta2 = np.zeros([interval,2])
x2 = np.linspace(0.,1.,interval)

for i in range(interval):
    beta2[i,0] = x2[i]
    beta2[i,1] = ((Sm * (1 - am)) / 4 / sigma / (1 - (1 - x2[i]) / 2)) ** (1./4.)
###END OF FILLING ARRAYS

#######################################################
### PLOT CONFIG

plt.plot(beta1[:,0], beta1[:,1], label = 'Earth')
plt.vlines(beta1[23,0], 0, beta1[23,1])
plt.plot(beta2[:,0], beta2[:,1], color = 'blue', linestyle = 'dashed', label = 'Mars')
plt.vlines(beta2[94,0], 0, beta2[94,1])
plt.legend(loc=1)
plt.title('"Surface Temperature as a Function of Beta"')
plt.xlabel('Beta')
plt.ylabel('Surface Temperature (K)')
plt.axis([0,1,200,320])
plt.xticks(np.arange(0, 1.1, .1))
### END OF PLOT CONFIG

plt.show()
