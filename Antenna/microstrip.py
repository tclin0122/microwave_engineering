import numpy as np
import cmath
import scipy.constants as sp
from scipy.special import jv

## parameter
# (m,n) = (1,1)
fc = 2.45 * np.power(10,9) # Hz
h = 1.52 * 0.001 # mm
dielectric = 3.38

Knm = 1.841

E0 = 1 # V/m
n = 1
wave_length = sp.speed_of_light / fc
k0 = 2 * np.pi / wave_length
k = k0 * np.sqrt (dielectric)

## Calculate
a = Knm  * sp.speed_of_light / (2 * np.pi * fc * np.sqrt(dielectric))
a_eff = a * np.sqrt (1 + (2 * h / (np.pi*a*dielectric)) * (np.log(np.pi*a/(2*h) + 1.7726)))
f_nm = Knm * sp.speed_of_light / (2 * np.pi * a_eff * np.sqrt(dielectric))

a_corr_eff = Knm * sp.speed_of_light / (2 * np.pi * fc * np.sqrt(dielectric))
a_g = 18.9 / 1000
a_corr =  a_g * np.sqrt (1 + (2 * h / (np.pi*a_g*dielectric)) * (np.log(np.pi*a_g/(2*h) + 1.7726)))
f_nm_corr = Knm * sp.speed_of_light / (2 * np.pi * a_corr * np.sqrt(dielectric))

print(a*1000)
print(a_eff*1000)
print(f_nm / np.power(10,9))

## Calculate Field
rho = 0
phi = 0
Ez = E0 * jv(n,k * rho) * np.cos(n*np.deg2rad(phi))
print(Ez)

rho = a
phi = 0
Ez = E0 * jv(n,k * rho) * np.cos(n*np.deg2rad(phi))
print(Ez)

## Far field
r = 10 #m
theta = 0
phi = 0

E_theta = ((np.power(1j,n)*h*a*k0*E0*jv(n,k*a)*cmath.exp(-1j*k0*r))/(2*r))* np.cos(n*np.deg2rad(phi)) * (jv(n+1,k0*a*np.sin(theta)) - jv(n-1,k0*a*np.sin(theta)))
E_phi = ((np.power(1j,n)*h*a*k0*E0*jv(n,k*a)*cmath.exp(-1j*k0*r))/(2*r))* np.cos(np.deg2rad(theta)) * np.sin(n*np.deg2rad(phi)) * (jv(n+1,k0*a*np.sin(theta)) + jv(n-1,k0*a*np.sin(theta)))

E = np.sqrt(E_theta*E_theta + E_phi*E_phi)
E = E_theta + E_phi
print (abs(E)*1000)