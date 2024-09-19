import numpy as np
import cmath
import scipy.constants as sp
from scipy.special import j1

## Parameters by the questions
Gt_dB = 35 #dBi
Pt = 100 # W
fc = 12 *np.power(10,9) #Hz
Pr_dBm = -85 #dBm
r = 35786 #km

## Trabsform the units
Pt = Pt * 1000 #mW
Pt_dBm = 10*np.log10(Pt)
r = r * 1000 #m
wave_length = sp.speed_of_light / fc #m
fspl_dB = 10 * np.log10(np.power(wave_length,2)/(np.power(4*np.pi*r,2))) 

## Link budget calculation
Gr_dB =  Pr_dBm - Pt_dBm - Gt_dB - fspl_dB
print("The minimum required antenna gain",Gr_dB, "dBi")

## Antenna calculation
a = (wave_length / (2*np.pi)) * np.power(10,(Gr_dB/20))
print("The minimum required diameter of the circular aperture.",2*a,"m")
HPBW = 1.02*wave_length/(2*a)
print("The beamwidth (in degrees) of the reflector antenna.",np.rad2deg(HPBW),"deg")
a = (wave_length / (2*np.pi)) * np.power(10,(Gr_dB/20)) / (np.sqrt(0.56))
print("The minimum required diameter of the circular aperture when p =2 is ",2*a,"m")
HPBW = 42.1*wave_length/(a)
print("The beamwidth (in degrees) of the reflector antenna.",HPBW,"deg")

# Next
print("")
print("second stage")
## Parameters
D = 4 #m
fc = 3 * np.power(10,9) #Hz
E0 = 100 # V/m
## Convert
wave_length = sp.speed_of_light / fc # m
a = D / 2 #m

## Calculate
HPBW = 29.2*wave_length/(a)
print("The beamwidth (in degrees) of the reflector antenna.",HPBW,"deg")
Rff = 2 * np.power(D,2) / wave_length
print("Fraunhofer distance where the far-field region starts at",Rff, "m")

## Calculate field
theta = HPBW / 2 # degree
phi = 45 # degree
r = 10000 # m

k0 = 2 * np.pi / wave_length
E_theta = ((1j*a*a*k0*E0*cmath.exp(-1j*k0*r))/(2*r))*(1+np.cos(np.deg2rad(theta)))*np.cos(np.deg2rad(phi))*(j1(k0*a*np.sin(np.deg2rad(theta))))/(k0*a*np.sin(np.deg2rad(theta)))
E_phi = ((-1j*a*a*k0*E0*cmath.exp(-1j*k0*r))/(2*r))*(1+np.cos(np.deg2rad(theta)))*np.sin(np.deg2rad(phi))*(j1(k0*a*np.sin(np.deg2rad(theta))))/(k0*a*np.sin(np.deg2rad(theta)))
print("Electric field strength of the θ component", abs(E_theta),"V/m ")
print("Electric field strength of the ϕ component", abs(E_phi),"V/m ")

# calculate power density usinging poynting vector
P = 0.5/ 377 * (E_theta * E_theta + E_phi * E_phi)
print("The power density (in W/m^2) at a distance r=10 km from the antenna with θ = HPBW/2 and ϕ=45 degrees. ", abs(P), "W/m^2")

## Calculate field
theta = HPBW / 1000000000 # degree
phi = 0 # degree
r = 10000 # m

k0 = 2 * np.pi / wave_length
E_theta = ((1j*a*a*k0*E0*cmath.exp(-1j*k0*r))/(2*r))*(1+np.cos(np.deg2rad(theta)))*np.cos(np.deg2rad(phi))*(j1(k0*a*np.sin(np.deg2rad(theta))))/(k0*a*np.sin(np.deg2rad(theta)))

E_phi = ((-1j*a*a*k0*E0*cmath.exp(-1j*k0*r))/(2*r))*(1+np.cos(np.deg2rad(theta)))*np.sin(np.deg2rad(phi))*(j1(k0*a*np.sin(np.deg2rad(theta))))/(k0*a*np.sin(np.deg2rad(theta)))

# calculate power density usinging poynting vector
P = 0.5/ 377 * (E_theta * E_theta + E_phi * E_phi)
print("The power density (in W/m^2) at a distance r=10 km from the antenna at broadside ", abs(P), "W/m^2")


# New section
## Parameter
D = 0.6 # m
a = D / 2
F = 0.5 * D # m
Zr = (D*D)/(16*F)
print("The depth of the reflector",Zr*100, "cm")