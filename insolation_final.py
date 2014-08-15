import math
import numpy as np
import matplotlib.pyplot as plt

#######################################################
### PARAMETERS

latitude1 = input("Please enter the latitude of the location #1:  \n")
latitude2 = input("Please enter the latitude of the location #2:  \n")
phi1 = math.pi * latitude1 / 180. # Location 1
phi2 = math.pi * latitude2 / 180. # Location 2
h = 0 # solar time in hours
sigma = 5.67e-8 # W / m^2 / K^4
L = 3.9e26 # W
solar_cons = L / (4. * math.pi * 1.5e11 ** 2) # W / m^2

### Fourier Coefficients

# Coefficients for lil_delta:

a1 = [0.006918, -0.399912, -0.006758, -0.002697]
b1 = [0, 0.070257, 0.000907, 0.001480]

# Coefficients for distance ratio:
a2 = [1.000110, 0.034221, 0.000719]
b2 = [0, 0.001280, 0.000077]
### END OF PARAMETERS

#######################################################
### INITIALIZED ARRAYS

theta_d = np.zeros(365)
lil_delta = np.zeros(365)
dist_ratio = np.zeros(365)

# Location 1 arrays
cos_zenith1 = np.zeros(365)
zenith1 = np.zeros(365)
insolation1 = np.zeros(365)

cos_zenith1_h = np.array(24)
zenith1_h = np.array(24)
insolation1_h = np.array(24)

# Location 2 arrays
cos_zenith2 = np.zeros(365)
zenith2 = np.zeros(365)
insolation2 = np.zeros(365)
### END OF INITIALIZED ARRAYS

#######################################################
### DAY LOOP
for day in range(365):
	theta_d[day] = 2 * math.pi * day / 365. # in radians

	# Fourier loop -- DECLANATION
	for n in range(len(a1)):
		lil_delta[day] += a1[n] * math.cos(n * theta_d[day]) + b1[n] * math.sin(n * theta_d[day])	 
	
	# Fourier loop -- DISTANCE RATIO
	for m in range(len(a2)):
		dist_ratio[day] += a2[m] * math.cos(m * theta_d[day]) + b2[m] * math.sin(m * theta_d[day])
	
	# Filling of Arrays
	cos_zenith1[day] = math.sin(phi1) * math.sin(lil_delta[day]) + math.cos(phi1) * math.cos(lil_delta[day]) * math.cos(h)
	zenith1[day] = math.acos(cos_zenith1[day]) * 180. / math.pi
	insolation1[day] = solar_cons * dist_ratio[day] ** 2. * cos_zenith1[day]
		
	cos_zenith2[day] = math.sin(phi2) * math.sin(lil_delta[day]) + math.cos(phi2) * math.cos(lil_delta[day]) * math.cos(h)
	zenith2[day] = math.acos(cos_zenith2[day]) * 180 / math.pi
	insolation2[day] = solar_cons * dist_ratio[day] ** 2 * cos_zenith2[day]
### END OF DAY LOOP

#######################################################
### INSOLATION ON DAY 172; JUNE 21 (SOLTICE)

day = 172

theta_d = 2 * math.pi * day / 365. # in radians
cos_zenith_h = cos_zenith1[day]

insolation1_h = np.zeros(24)

# Fourier loop -- DECLANATION
lil_delta_h = 0
for n in range(len(a1)):
	lil_delta_h += a1[n] * math.cos(n * theta_d) + b1[n] * math.sin(n * theta_d)	 

# Fourier loop -- DISTANCE RATIO
dist_ratio_h = 0
for m in range(len(a2)):
	dist_ratio_h += a2[m] * math.cos(m * theta_d) + b2[m] * math.sin(m * theta_d)

for hour in range(24):
    h = 2 * hour * math.pi / 24 
    cos_zenith1_h = math.sin(phi1) * math.sin(lil_delta_h) + math.cos(phi1) * math.cos(lil_delta_h) * math.cos(h)
    zenith1_h = math.acos(cos_zenith1_h) * 180 / math.pi
    insolation1_h[hour] = solar_cons * dist_ratio_h ** 2 * cos_zenith1_h

### END OF INSOLATION ON DAY 172; JUNE 21 (SOLTICE)

#######################################################
### OUTPUT TO INSOLATION1.TXT

# Create output file
out = open('insolation1.txt', 'w')

# Convert insolation1 to string and write to file
outputstr = ['\n']*len(insolation1)
for i in range(len(insolation1)):
    outputstr[i] = str(insolation1[i]) + outputstr[i]
out.writelines(outputstr)
out.close()
print outputstr
### END OF OUTPUT TO INSOLATION1.TXT
#------------------------------------------------------

### OUTPUT TO INSOLATION2.TXT

# Create output file
out = open('insolation2.txt', 'w')

# Convert insolation1 to string and write to file
outputstr = ['\n']*len(insolation2)
for i in range(len(insolation2)):
    outputstr[i] = str(insolation2[i]) + outputstr[i]
out.writelines(outputstr)
out.close()
print outputstr
### END OF OUTPUT TO INSOLATION1.TXT

#######################################################
### PLOT CONFIG 

plt.figure(1)
plt.plot(insolation1, label = '34 deg. N')
plt.plot(insolation2, color='blue', linestyle = 'dashed', label = '3 deg. S')
plt.legend(loc=8)
plt.title('Daily Insolation at Two Latitudes')
plt.xlabel('Day of Year')
plt.ylabel('Insolation (W/m^2)')
plt.axis([0,364,700,1500])
### END OF PLOT CONFIG
#------------------------------------------------------

plt.figure(2)
plt.title('Insolation at 34 deg. Latitude, on June 21')
plt.scatter(range(24),insolation1_h, c = 'r', s=100)
plt.xlabel('Solar Time (hrs.)')
plt.axis([0,23,0,1300])
plt.ylabel('Insolation (W/m^2)')

plt.show()
### END PLOT CONFIG

#######################################################
"""
### DISCUSSION

# Part II:
	The plots of insolation at 34 deg. N and 3 deg. S show that the intensity of
the Sun's radiation is different at different points on Earth at different times
of the year. This is expected because there is an inverse square relation between
the Sun's flux density of radiated energy and the distance to some point on Earth.
Since Earth has spherical symmetry, there is a shared zenith angle for all points
along a given latitude line so there is no need to regard longitude.
	As the day of the year varies, the Earth moves through its elliptic orbit of 
the Sun. The insolation plots vs. time of year indicate the amount of available
radiation during a specific season. In the Northern Hemisphere, maximum insolation
coincides with movement throught the apehelion adn occurs during the months of June, 
July, and August. Within the tropics, locations along the 3 deg. latitude receive, 
on average, greater amounts of radiation since the Sun can have a larger declanation.
This is why there is far more radioation reaching the tropical location on the
equinox. 

Part III:
	Insolation for a given latitude varies by hours of the day. This graph shows that 
models the rising and setting of the Sun on the Summer Solstice.
"""
