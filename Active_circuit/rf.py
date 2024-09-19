import numpy as np
import math
import scipy.constants as sy

## S parameters
Mag_S11 = 0.65
Ang_S11 = -93
Mag_S12 = 0.04
Ang_S12 = 47
Mag_S21 = 17.88
Ang_S21 = 116
Mag_S22 = 0.65
Ang_S22 = -46


## Convert S parameters to complex numbers
S11 = complex(Mag_S11 * np.cos(np.deg2rad(Ang_S11)),Mag_S11 * np.sin(np.deg2rad(Ang_S11)))
S21 = complex(Mag_S21 * np.cos(np.deg2rad(Ang_S21)),Mag_S21 * np.sin(np.deg2rad(Ang_S21)))
S12 = complex(Mag_S12 * np.cos(np.deg2rad(Ang_S12)),Mag_S12 * np.sin(np.deg2rad(Ang_S12)))
S22 = complex(Mag_S22 * np.cos(np.deg2rad(Ang_S22)),Mag_S22 * np.sin(np.deg2rad(Ang_S22)))

Z0 =50 #impedance
Zload = 50.000000001
Zsource = 100
ZL = (Zload-Z0) / (Zload+Z0)
ZS = (Zsource-Z0) / (Zsource+Z0)


refl_input = S11 + (S12*S21)/(1/ZL - S22)
print("The absolute value of the input reflection coefficient = ",abs(refl_input))

refl_output = S22 + (S12*S21)/(1/ZS - S11)
print("The absolute value of the output reflection coefficient = ",abs(refl_output))


## Gain calculation
GT = ((1 - np.power(abs(ZL),2))/(np.power(abs(1 - S22*ZL),2))) * np.power(abs(S21),2) * ((1 - np.power(abs(ZS),2))/(np.power(abs(1 - refl_input*ZS),2)))
GT_dB = 10*np.log10(GT)
Gav = ((1 - np.power(abs(ZS),2))/(np.power(abs(1 - S11*ZS),2))) * np.power(abs(S21),2) * (1/(1 - np.power(abs(refl_output),2)))
Gav_dB = 10*np.log10(Gav)
Gdel = ((1 - np.power(abs(ZL),2))/(np.power(abs(1 - S22*ZL),2))) * np.power(abs(S21),2) * (1/(1 - np.power(abs(refl_input),2)))
Gdel_dB = 10*np.log10(Gdel)

print("Transducer power gain",GT_dB,"dB")       # mismatch at input and output
print("Available power Gain",Gav_dB,"dB")       # mismatch at input
print("Delivered power Gain ",Gdel_dB,"dB")     # mismatch at output

## Noise power at the input of the amplifier
B = 9.548 * np.power(10,6)
k = sy.Boltzmann
T = 290
Pn = k*T*B
print("Pn",Pn)
print(k)

Vs = 0.1/1000 *1/np.sqrt(2)
Rs = 100
Pn_dB = 10*np.log10(Pn*1000)
print(Pn_dB)

Ps = np.power(Vs,2)/(4*Rs)
Ps_dB = 10*np.log10(Ps/0.001)
print("Available power from the source",Ps_dB,"dBm")
print("Available power at the output of the amplifier",Ps_dB + GT_dB,"dBm")
SNR = Ps_dB + GT_dB - Pn_dB 
print(SNR)

