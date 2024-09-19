import numpy as np
import cmath
import scipy.constants as sp
from scipy.special import jv

## parameter
# (m,n) = (1,1)
fc = 5 * np.power(10,9) # Hz
h = 1.52 * 0.001 # mm
dielectric = 3.38
wave_length = sp.speed_of_light / fc
k0 = 2 * np.pi / wave_length
k = k0 * np.sqrt (dielectric)
print(k0)

Knm = 1.841
ak = 1.0 # uniform rapering
K = 16 # 16 element in the array

## Calculate
d = wave_length / 2
print("The element spacing ",d*1000, "mm")

a = Knm  * sp.speed_of_light / (2 * np.pi * fc * np.sqrt(dielectric))

print("The radius of the microstrip antenna. (Neglect fringing-field effects) ",a*1000, "mm")

## Radiation pattern
theta = 45
scan_ang = 0
i = 1
s_tot = 0.0
while (i < K+1):
    s_tot = s_tot + ak * np.exp(1j*k0*(i-1)*d*(np.sin(np.deg2rad(theta) - np.sin(np.deg2rad(scan_ang)))))
    i = i + 1
    
            
# S = (jv(2,k0*a*np.sin(np.deg2rad(theta))) - jv(0,k0*a*np.sin(np.deg2rad(theta)))) * s_tot
S = (jv(2,k0*a*np.sin(np.deg2rad(theta))) - jv(0,k0*a*np.sin(np.deg2rad(theta))))

print(20*np.log10(abs(S)))