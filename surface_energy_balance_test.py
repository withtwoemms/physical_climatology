"""
This program extracts gradient data compiled by a sensor array in Williamsport, PA and calculates incoming and outgoing long-wave radiation, latent and sensible heat fluxes, ground heat flux and net radiation at the surface. Wind profiles are generate from the same data and the 'roughness length' is calculated. 
"""

import numpy as np
from numpy import genfromtxt
import math
import matplotlib.pyplot as plt

#######################################################
### GLOBAL CONSTANTS

Cp = 1004 # air heat capacity, J/kg/K
k = 0.4 # von Karmen
kb = 5.67 # Stefan-Boltzmann Constant, W/m^2
z = 2 # sensor height in m
rho_a = 1 # air density, kg/m^3
g = 9.81 # m/s^2 
epsilon = 0.95 # emissivity
sigma = 5.67 * 10**(-8) # Stefan-Boltzmann constant, W/m^2/K^4

#######################################################
### DATA IMPORT

input_data = genfromtxt('lycoming2001082_f03.txt', delimiter='\t')

hours = input_data[:,2]
G = input_data[:,15]
SW_DOWN = input_data[:,13] 
SW_UP = input_data[:,14]

### data @ 3m
temp_3m_celsius = input_data[:,4]
temp_3m_kelvin = temp_3m_celsius + 273
rh_3m_percent = input_data[:,5]
wind_speed_3m = input_data[:,6]
wind_dir_3m = input_data[:,7] 

### data @ 1m
temp_1m_celsius = input_data[:,9]
temp_1m_kelvin = temp_1m_celsius + 273
rh_1m_percent = input_data[:,10]
wind_speed_1m = input_data[:,11]
wind_dir_1m = input_data[:,12]

#######################################################
### CALCULATION OF H, LE

### initialization of arrays
es_3m = np.zeros(len(hours)) # saturation vapor pressure in Pascals
es_1m = np.zeros(len(hours))
e_3m = np.zeros(len(hours)) # vapor pressure in Pascals
e_1m = np.zeros(len(hours)) 
rho_v_3m = np.zeros(len(hours)) # vapor density @ 3m in kg/m^3
rho_v_1m = np.zeros(len(hours)) # vapor density @ 3m in kg/m^3
dt_dz = np.zeros(len(hours)) # temperature gradient K/m
drho_dz = np.zeros(len(hours)) # water vapor density gradient kg/m^4 
temp_kelvin_avg = np.zeros(len(hours)) # average temperature over dz
Ri = np.zeros(len(hours)) # Richardson Number
phi_m = np.zeros(len(hours)) 
K_ratio = np.zeros(len(hours))
du_dz = np.zeros(len(hours)) # vertical wind speed profile in Hz
H = np.zeros(len(hours)) # sensible heat
LE = np.zeros(len(hours)) # latent heat
Q_net = np.zeros(len(hours)) # net radiation

tau = np.zeros(len(hours)) # shearing stress
u_star = np.zeros(len(hours)) # friction velocity
z_naught = np.zeros(len(hours)) # roughness length

### calculations
for i in range (len(hours)):
    es_3m[i] = 611.2 * np.exp(17.67 * temp_3m_celsius[i] / (temp_3m_celsius[i] + 243.5))
    es_1m[i] = 611.2 * np.exp(17.67 * temp_1m_celsius[i] / (temp_1m_celsius[i] + 243.5))
    e_3m[i] = rh_3m_percent[i] / 100 * es_3m[i]
    e_1m[i] = rh_1m_percent[i] / 100 * es_1m[i]
    rho_v_3m[i] = 2.17 * e_3m[i] / 1000 / temp_3m_kelvin[i]
    rho_v_1m[i] = 2.17 * e_1m[i] / 1000 / temp_1m_kelvin[i]
    drho_dz[i] = (rho_v_3m[i] - rho_v_1m[i]) / z    

    dt_dz[i] = (temp_3m_kelvin[i] - temp_1m_kelvin[i]) / z
        
    temp_kelvin_avg[i] = (temp_3m_kelvin[i] + temp_1m_kelvin[i]) / z

    Ri[i] = g / temp_kelvin_avg[i] * (((temp_3m_kelvin[i] - temp_1m_kelvin[i]) / z) / ((wind_speed_3m[i] - wind_speed_1m[i]) / z)**2)

    du_dz[i] = (wind_speed_3m[i] - wind_speed_1m[i]) / z

    if Ri[i] < 0:
        phi_m[i] = (1 - 16 * Ri[i])**(-1 / 3.)
        K_ratio[i] = (1 - 16 * Ri[i])**(1 / 4.)
    elif Ri[i] == 0:
        phi_m[i] = 1
        K_ratio[i] = 1
    elif Ri[i] > 0:
        phi_m[i] = (1 + Ri[i])**(1 / 3.)
        K_ratio[i] = 1
    
    H[i] = -1 * rho_a * Cp * dt_dz[i] * K_ratio[i] * du_dz[i] * k**2 * z**2 * phi_m[i]**(-2)
    
    tau[i] = rho_a * (du_dz[i] * k * z / phi_m[i])**2
    u_star[i] = (tau[i] / rho_a)**(1 / 2.)
    z_naught[i] = math.exp(-(k * wind_speed_1m[i]) / (u_star[i] * phi_m[i]))


    if temp_kelvin_avg[i] > 273:
       L = 2.5 * 10**6 # J / kg
    elif temp_kelvin_avg[i] < 273:
        L = 2.83 * 10**6 
    LE[i] = -1 * L * drho_dz[i] * K_ratio[i] * du_dz[i] * k**2 * z**2 * phi_m[i]**(-2)

#######################################################
### CALCULATION OF LW_DOWN, LW_UP, Q*, G

LW_DOWN = np.zeros(len(hours)) # incoming longwave
LW_UP = np.zeros(len(hours)) # outgoing longwave

for i in range(len(hours)):
    LW_DOWN[i] = 5.31 * 10**(-13) * temp_kelvin_avg[i]**6
    LW_UP[i] = (1 - epsilon) * LW_DOWN[i] + epsilon * sigma + temp_1m_kelvin[i]

    Q_net[i] = (LW_DOWN[i] - LW_UP[i]) + (SW_DOWN[i] - SW_UP[i])
    
    G[i] = Q_net[i] - H[i] - LE[i]

#######################################################
### CALCULATION OF SURFACE ALBEDO

alpha = np.zeros(len(hours))

for i in range(len(hours)):
    alpha[i] = 1 - (SW_DOWN[i] - SW_UP[i]) / SW_DOWN[i]

"""
Some values return 'nan' due to division by zero. This is a non-physical answer most likely due to insensitvity of the incoming solar radiation sensor
"""
#######################################################
### CALCULATION OF WIND PROFILES; ROUGHNESS LENGTH

U6 = [0, wind_speed_1m[6], wind_speed_3m[6]]
U13 =  [0, wind_speed_1m[13], wind_speed_3m[13]]

ln_U6 = [wind_speed_1m[6], wind_speed_3m[6]]
ln_U13 = [wind_speed_1m[13], wind_speed_3m[13]]

### is the atmosphere stable?

if phi_m[6] < 0:
    stability_6am = 'unstable at 6 am'
elif phi_m[6] > 0: 
    stability_6am = 'stable at 6 am'
elif phi_m[6] == 0:
    stability_6am = 'neutral at 6 am'

if phi_m[13] < 0:
    stability_6am = 'unstable at 1 pm'
elif phi_m[13] > 0: 
    stability_1pm = 'stable at 1 pm'
elif phi_m[13] == 0:
    stability_1pm = 'neutral at 1 pm'

print stability_6am, stability_1pm

#######################################################
### PLOT CONFIGS

x1 = np.zeros(len(hours))
for i in range(len(hours)):
    x1[i] = i + 1

x2 = [0., 1., 3.]
ln_x2 = [math.log(x2[1]), math.log(x2[2])]

### PLOT OF Q, LW_UP, LW_DOWN, SW_UP, SW_DOWN
plt.figure(1)
plt.title('Hourly Energy Flux: Williamsport, PA; Day 82, 2001')
plt.xlabel('Hour')
plt.ylabel('Energy Flux (W/m^2)')
plt.plot(x1, Q_net, label = 'Net Radiation')
plt.plot(x1, LW_UP, ls = '-', label = 'outgoing LW')
plt.plot(x1, LW_DOWN, ls = '--', label = 'incoming LW')
plt.plot(x1, SW_UP, ls = '-.', label = 'outgoing solar')
plt.plot(x1, SW_DOWN, ls = ':', label = 'incoming solar')
plt.legend(loc=2)
plt.axis([1,24,-100,750])

### PLOT OF Q, H, LE, G, ALBEDO
plt.figure(2)
plt.title('Energy Balance: Williamsport, PA; Day 82, 2001')
plt.xlabel('Hour')
plt.ylabel('Energy Flux (W/m^2)')
plt.plot(x1, Q_net, label = 'Net Radiation')
plt.plot(x1, H, ls = '--', label = 'Sensible Heat')
plt.plot(x1, LE, ls = ':', label = 'Latent Heat')
plt.plot(x1, G, ls = '-.', label = 'Ground Heat')
plt.legend(loc=2)

### surface albedo inset
a = plt.axes([.2, .3, .2, .2])
plt.scatter(x1, alpha)
plt.title('Surface Albedo')
plt.setp(a, xticks=[5,10,15,20], yticks=[0.05,0.1,0.15,.2,.25,.3,.35])

#plt.axis([1,24,-2150,500])

### PLOT OF WIND PROFILE
plt.figure(3)
plt.title('Wind Profile: Williamsport, PA; Day 82, 2001')
plt.xlabel('Wind Speed, U(z)')
plt.ylabel('Height, z')
plt.plot(U6, x2, label = '6 am')
plt.plot(U13, x2, ls = '-.', label = '1 pm')
plt.legend(loc=4)
plt.text(2.75, 1.5, r'$U(z) = \left ( \frac{u_\star \phi_m}{k}  \right ) ln\left (\frac{z_0}{z} \right )$', fontsize = 15)

### LOG PLOT OF WIND PROFILE + FITS
fig4 = plt.figure(4)
plt.title('Wind Profile: Williamsport, PA; Day 82, 2001')
plt.xlabel('Wind Speed, U(z)')
plt.ylabel('Height, ln(z)')
plt.text(-0.8, 2.5, r'$\left (\frac{k}{u_\star  \phi_m}  \right ) U(z) + ln(z_0) = ln(z)$', fontsize = 15)
plt.scatter(ln_U6, ln_x2, label = '6 am', s = 60)
plt.scatter(ln_U13, ln_x2, label = '1 pm', marker = '+', s = 80)

### linear fits
xfit1 = np.linspace(0, 3, 30)
xfit2 = np.linspace(0, 7, 70) 
fit1 = np.polyfit(ln_U6, ln_x2, 1)
fit2 = np.polyfit(ln_U13, ln_x2, 1)
fit_fn1 = np.poly1d(fit1)
fit_fn2 = np.poly1d(fit2)

plt.plot(xfit1, fit_fn1(xfit1),'-.')
plt.plot(xfit2, fit_fn2(xfit2), ':')
plt.legend(loc=1)

### zoom inset
b = plt.axes([.6, .15, .3, .3])
plt.plot(xfit1[0:2], fit_fn1(xfit1)[0:2],'-.')
plt.plot(xfit2[0:2], fit_fn2(xfit2)[0:2], ':')
plt.title('Zoom')
plt.setp(b, xticks=[0,0.05,.1], yticks=[-3.25,-3.2,-3.15,-3.1,-3.05])

plt.show()

#######################################################
### DISCUSSION
"""
### Figure 1: Net, LW, and SW Radiation Fluxes
    There is a net positive flux of energy from the surface. The flux reaches a maximum at 1 pm, drops dramatically by 2 pm, and recovers significantly by 3 pm. This perculiarity is the most remarkable feature of the plot and may be explained by the possible passing of a cloud. Such a hypothesis is reasonable because insolation is similarly reduced. The same pattern is borne by the sensible and latent heat flux curves of figure 2. Interestingly, the latent heat curve remains relatively constant which would be in line with a pause in evaporation. Such a pause would be expected in a scenario where there was a period of reduced insolation such as in the event of a cloud passage.

### Figure 2: Net, Sensible, Latent and Ground Energy Fluxes
    The data suggest that Williamstown, PA is relatively dry with a cold surface during the month of March. The sensible heat flux (H) correlates very well with the pattern outgoing solar radiation (see figure 1). Latent heat (LE) is the lowest flux with a maximum of ~40 W/m^2 at noon with the expected negative fluxes in the morning and evening most likely attributable to the formation of dew. The low latent heat flux suggests a dry surface. From the ground heat flux (G), it can be inferred that the surface is cold, i.e. that there is a strong, negative temperature gradient between the surface and sub-surface, since the ground heat flux is nearly double that of sensible heat at noon. There exists an anomolous dip in the ground heat flux at peak insolation most likely due to reduction of the the surface/sub-surface temperature gradient.  

### Figure 3: Wind Profile
    Wind speeds at 6 am approximately triple by 1 pm. Increased convection due to conductive heat exchange between the surface and atmosphere is a possible reason for why the wind speed profiles differ. By 1 pm, surface insolation is greater than at 6 am.
    Accurate conclusions about the stability of the atmosphere cannot be drawn with only two data points. Having wind speed measurements at higher altitudes can show whether or not the natural log of the profile significantly deviates from linearity and thus would make it possible to diagnose the stability. The positive values for the phi parameter indicate that the atmosphere is stable. A stable atmosphere has horizontally, oblong-shaped eddy currents which are most efficient in transporting latent heat. 

### Figure 4: Wind Profile
    The 'roughness lengths' are the exponentials of the y-intercept values featured in figure 4. Values of approximately 3.88 and 4.16 cm were determined for wind speed profiles taken at 6 am and 1 pm, respectively. These valuse are consistent with a generally flat surface most likely covered with a grass. The presence of some flexible foliage such as grass explains why the roughness length increases marginally later in the day when the wind speed profile spans greater speeds. A lower wind speed requires less mass to oppose do winds of higher speed. 
"""
