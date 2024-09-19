import numpy as np
import math
import scipy.constants as sy

## S parameters
Mag_S11 = 0.47
Ang_S11 = 173
Mag_S12 = 0.06
Ang_S12 = 30
Mag_S21 = 6.01
Ang_S21 = 55
Mag_S22 = 0.27
Ang_S22 = -75

## Convert S parameters to complex numbers
S11 = complex(Mag_S11 * np.cos(np.deg2rad(Ang_S11)),Mag_S11 * np.sin(np.deg2rad(Ang_S11)))
S21 = complex(Mag_S21 * np.cos(np.deg2rad(Ang_S21)),Mag_S21 * np.sin(np.deg2rad(Ang_S21)))
S12 = complex(Mag_S12 * np.cos(np.deg2rad(Ang_S12)),Mag_S12 * np.sin(np.deg2rad(Ang_S12)))
S22 = complex(Mag_S22 * np.cos(np.deg2rad(Ang_S22)),Mag_S22 * np.sin(np.deg2rad(Ang_S22)))

## calculate stability
Zload = 50.000000001
Zsource = 50.000000001
Z0 =50
ZL = (Zload-Z0) / (Zload+Z0)
ZS = (Zsource-Z0) / (Zsource+Z0)
refl_input = S11 + (S12*S21)/(1/ZL - S22)
refl_output = S22 + (S12*S21)/(1/ZS - S11)
#print(abs(refl_input))
#print(abs(refl_output))
Delta = abs(S11*S22 - S21*S12)
K = (1 - np.power(abs(S11),2) - np.power(abs(S22),2) + np.power(abs(Delta),2))/(2*abs(S21*S12))
mu = (1 - np.power(abs(S22),2))/(abs(S11-np.conj(S22)*Delta) + abs(S12*S21))
print("K = ",K)
print("∆ = ",Delta)
print("µ =",mu)

## Calculate power gain
Gav = ((1 - np.power(abs(ZS),2))/(np.power(abs(1 - S11*ZS),2))) * np.power(abs(S21),2) * (1/(1 - np.power(abs(refl_output),2)))
Gav_dB = 10*np.log10(Gav)

print("Available power gain", Gav_dB, "dB")


## S parameters
Mag_S11 = 0.869
Ang_S11 = -159
Mag_S12 = 0.031
Ang_S12 = -9
Mag_S21 = 4.250
Ang_S21 = 61
Mag_S22 = 0.507
Ang_S22 = -117

## Convert S parameters to complex numbers
S11 = complex(Mag_S11 * np.cos(np.deg2rad(Ang_S11)),Mag_S11 * np.sin(np.deg2rad(Ang_S11)))
S21 = complex(Mag_S21 * np.cos(np.deg2rad(Ang_S21)),Mag_S21 * np.sin(np.deg2rad(Ang_S21)))
S12 = complex(Mag_S12 * np.cos(np.deg2rad(Ang_S12)),Mag_S12 * np.sin(np.deg2rad(Ang_S12)))
S22 = complex(Mag_S22 * np.cos(np.deg2rad(Ang_S22)),Mag_S22 * np.sin(np.deg2rad(Ang_S22)))

## calculate stability
Zload = 50.000000001
Zsource = 50.000000001
Z0 =50
ZL = (Zload-Z0) / (Zload+Z0)
ZS = (Zsource-Z0) / (Zsource+Z0)
refl_input = S11 + (S12*S21)/(1/ZL - S22)
refl_output = S22 + (S12*S21)/(1/ZS - S11)
#print(abs(refl_input))
#print(abs(refl_output))
Delta = abs(S11*S22 - S21*S12)
K = (1 - np.power(abs(S11),2) - np.power(abs(S22),2) + np.power(abs(Delta),2))/(2*abs(S21*S12))
mu = (1 - np.power(abs(S22),2))/(abs(S11-np.conj(S22)*Delta) + abs(S12*S21))
print("K = ",K)
print("∆ = ",Delta)
print("µ =",mu)

rs = abs(S12*S21/(np.power(abs(S11),2)-np.power(abs(Delta),2)))
cs = np.conj(S11 - Delta*np.conj(S22))/(np.power(abs(S11),2)-np.power(abs(Delta),2))

rl = abs(S12*S21/(np.power(abs(S22),2)-np.power(abs(Delta),2)))
cl = np.conj(S22 - Delta*np.conj(S11))/(np.power(abs(S22),2)-np.power(abs(Delta),2))

print("The radius of the output (load) stability circle = ", rl)
print("The centre of the output (load) stability circle.(real part) of the centre = ",cl, "wrong")
print("The radius of the input (source) stability circle. = ", rs)
print("The centre of the output (source) stability circle.(imaginary part) of the centre = ",cs, "wrong")
